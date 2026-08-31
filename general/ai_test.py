#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Interactive test Flask app that generates a sine-wave line plot.

Serves a single page in the browser with sliders for amplitude, frequency
and phase. Each slider change requests a freshly rendered plot image from
the server, to verify that a Flask + matplotlib environment is set up
correctly.
"""

# Standard library imports
import io
import string

# Third-party imports
import flask
import matplotlib
# force the non-interactive Agg backend before importing pyplot, since
# this process has no display and matplotlib must stay thread-safe for
# concurrent Flask requests
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Define variables
# =============================================================================
# host / port the development server listens on
HOST = '127.0.0.1'
PORT = 5000
# number of sample points across the x-axis range
NUM_POINTS = 400
# x-axis range, in radians
X_MIN, X_MAX = 0, 4 * np.pi
# single line color (muted blue, readable on a light surface)
LINE_COLOR = '#3B82C4'
# initial / min / max / step values for each sine wave parameter
# (initial, min, max, step)
AMPLITUDE_RANGE = (1.0, 0.1, 5.0, 0.1)
FREQUENCY_RANGE = (1.0, 0.1, 5.0, 0.1)
PHASE_RANGE = (0.0, -np.pi, np.pi, 0.01)

# HTML page: three range sliders and a plot image, refreshed via fetch()
# whenever a slider moves
INDEX_HTML = """
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Interactive Sine Wave</title>
  <style>
    body {
      font-family: -apple-system, sans-serif;
      background: #F7F8FA;
      color: #1F2937;
      display: flex;
      flex-direction: column;
      align-items: center;
      padding: 24px;
    }
    h1 { font-size: 18px; font-weight: 600; }
    #plot { max-width: 640px; width: 100%; }
    .controls {
      width: 100%;
      max-width: 480px;
      margin-top: 12px;
    }
    .control-row {
      display: flex;
      align-items: center;
      gap: 12px;
      margin: 8px 0;
    }
    .control-row label {
      width: 90px;
      font-size: 13px;
      color: #4B5563;
    }
    .control-row input[type="range"] {
      flex: 1;
      accent-color: $line_color;
    }
    .control-row span {
      width: 48px;
      text-align: right;
      font-size: 13px;
      color: #4B5563;
    }
  </style>
</head>
<body>
  <h1>Interactive Sine Wave</h1>
  <img id="plot" src="/plot.png?amplitude=$amp_init&frequency=$freq_init
&phase=$phase_init" alt="sine wave plot">
  <div class="controls">
    <div class="control-row">
      <label for="amplitude">Amplitude</label>
      <input type="range" id="amplitude" min="$amp_min" max="$amp_max"
             step="$amp_step" value="$amp_init">
      <span id="amplitude-value">$amp_init</span>
    </div>
    <div class="control-row">
      <label for="frequency">Frequency</label>
      <input type="range" id="frequency" min="$freq_min" max="$freq_max"
             step="$freq_step" value="$freq_init">
      <span id="frequency-value">$freq_init</span>
    </div>
    <div class="control-row">
      <label for="phase">Phase</label>
      <input type="range" id="phase" min="$phase_min" max="$phase_max"
             step="$phase_step" value="$phase_init">
      <span id="phase-value">$phase_init</span>
    </div>
  </div>
  <script>
    // re-request the plot image whenever any slider changes, reading
    // the current value of all three sliders each time
    const ids = ['amplitude', 'frequency', 'phase'];
    const sliders = Object.fromEntries(
      ids.map((id) => [id, document.getElementById(id)])
    );
    const valueLabels = Object.fromEntries(
      ids.map((id) => [id, document.getElementById(id + '-value')])
    );
    const img = document.getElementById('plot');

    function updatePlot() {
      const params = ids
        .map((id) => id + '=' + sliders[id].value)
        .join('&');
      img.src = '/plot.png?' + params;
      ids.forEach((id) => {
        valueLabels[id].textContent = sliders[id].value;
      });
    }

    ids.forEach((id) => {
      sliders[id].addEventListener('input', updatePlot);
    });
  </script>
</body>
</html>
"""

app = flask.Flask(__name__)


# =============================================================================
# Define functions
# =============================================================================
def sine_wave(x, amplitude, frequency, phase):
    """
    Compute a sine wave over x for the given parameters.

    :param x: array of x values, in radians
    :param amplitude: peak height of the wave
    :param frequency: number of oscillations per 2*pi of x
    :param phase: horizontal shift, in radians

    :return: array of y values, y = amplitude * sin(frequency * x + phase)
    """
    return amplitude * np.sin(frequency * x + phase)


def render_plot_png(amplitude, frequency, phase):
    """
    Render the sine wave for the given parameters to a PNG image.

    :param amplitude: peak height of the wave
    :param frequency: number of oscillations per 2*pi of x
    :param phase: horizontal shift, in radians

    :return: raw PNG image bytes
    """
    x = np.linspace(X_MIN, X_MAX, NUM_POINTS)
    y = sine_wave(x, amplitude, frequency, phase)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x, y, color=LINE_COLOR, linewidth=2)
    ax.set_xlabel('x (radians)')
    ax.set_ylabel('amplitude * sin(frequency * x + phase)')
    ax.set_ylim(-5.5, 5.5)
    # keep only left/bottom spines for a cleaner look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    # write directly to an in-memory buffer, no temp file needed
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png', dpi=110)
    # release the figure so repeated requests don't leak memory
    plt.close(fig)
    buffer.seek(0)
    return buffer.read()


# =============================================================================
# Define worker functions
# =============================================================================
def get_float_arg(name, default):
    """
    Read a float query-string argument from the current Flask request.

    :param name: query-string parameter name to read
    :param default: value to fall back to if the parameter is absent or
                     cannot be parsed as a float

    :return: parsed float value, or default
    """
    try:
        return float(flask.request.args.get(name, default))
    except (TypeError, ValueError):
        return default


@app.route('/')
def index():
    """
    Serve the main page with the plot image and its control sliders.

    :return: rendered HTML page as a string
    """
    # string.Template ($name) is used instead of str.format() because
    # the page's CSS/JS already contains many literal { } characters
    return string.Template(INDEX_HTML).substitute(
        line_color=LINE_COLOR,
        amp_init=AMPLITUDE_RANGE[0], amp_min=AMPLITUDE_RANGE[1],
        amp_max=AMPLITUDE_RANGE[2], amp_step=AMPLITUDE_RANGE[3],
        freq_init=FREQUENCY_RANGE[0], freq_min=FREQUENCY_RANGE[1],
        freq_max=FREQUENCY_RANGE[2], freq_step=FREQUENCY_RANGE[3],
        phase_init=PHASE_RANGE[0], phase_min=PHASE_RANGE[1],
        phase_max=PHASE_RANGE[2], phase_step=PHASE_RANGE[3],
    )


@app.route('/plot.png')
def plot_png():
    """
    Serve the sine wave plot as a PNG, using slider values from the
    request's query string.

    :return: Flask response with the PNG image bytes and image/png
             mimetype
    """
    amplitude = get_float_arg('amplitude', AMPLITUDE_RANGE[0])
    frequency = get_float_arg('frequency', FREQUENCY_RANGE[0])
    phase = get_float_arg('phase', PHASE_RANGE[0])
    png_bytes = render_plot_png(amplitude, frequency, phase)
    return flask.Response(png_bytes, mimetype='image/png')


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
    app.run(host=HOST, port=PORT, debug=True)

# =============================================================================
# End of code
# =============================================================================
