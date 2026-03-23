import re
from datetime import datetime
from datetime import timedelta

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy.time import Time
from tqdm import tqdm


# =============================================================================
# Define variables
# =============================================================================
# zoom in size in hours
ZOOM_IN_SIZE = 1.0
# define the recipe to look for in the log file
RECIPE = 'apero_preprocess_spirou'
# Define the regular expression for the line
PATTERN = (r"(?:\x1b\[[0-9;]*m)?(\d{2}:\d{2}:\d{2}\.\d{3})."
           r"*Recipe\s+([\w_]+)\s+has\s+been\s+successfully"
           r"\s+completed\s*[\t\s]+\((\d+\.\d+)\s+seconds\)")


# =============================================================================
# Define functions
# =============================================================================
def extract_log_to_table(file_path, start_date: str, recipe=None):
    # Updated Regex Breakdown:
    # (?:\x1b\[[0-9;]*m)? -> Non-capturing group to ignore ANSI color codes
    # (\d{2}:\d{2}:\d{2}\.\d{3}) -> Group 1: The timestamp
    # .*Recipe\s+([\w_]+) -> Group 2: The recipe name (optional, but useful!)
    # .*successfully\s+completed\s* -> Anchor text
    # \t?\((\d+\.\d+)\s+seconds\) -> Group 3: The duration

    pattern = re.compile(PATTERN)

    rows = []
    try:
        current_date = datetime.fromisoformat(start_date).date()
    except ValueError as exc:
        raise ValueError('start_date must be ISO format YYYY-MM-DD') from exc

    previous_time = None
    rollover_threshold_seconds = 12 * 3600

    # open file
    print(f'\nReading file {file_path}')
    with open(file_path, 'r') as file:
        lines = file.readlines()

    print(f'\nAnalysing {len(lines)} lines')
    # loop around lines
    for line in tqdm(lines):
        match = pattern.search(line)
        if match:
            ts_str, recipe_name, duration = match.groups()

            # Mask out recipes we don't want before updating date state.
            if recipe is not None and recipe_name != recipe:
                continue

            current_time = datetime.strptime(ts_str, '%H:%M:%S.%f').time()

            # Only treat large backward jumps as midnight rollover.
            if previous_time is not None and current_time < previous_time:
                prev_sec = (previous_time.hour * 3600 + previous_time.minute * 60
                            + previous_time.second
                            + previous_time.microsecond / 1e6)
                curr_sec = (current_time.hour * 3600 + current_time.minute * 60
                            + current_time.second
                            + current_time.microsecond / 1e6)
                if (prev_sec - curr_sec) > rollover_threshold_seconds:
                    current_date += timedelta(days=1)
            previous_time = current_time

            # Build a full timestamp using the rolling date anchor.
            full_ts = f"{current_date.isoformat()} {ts_str}"


            # append to rows
            rows.append({
                'timestamp': full_ts,
                'recipe': recipe_name,
                'runtime': float(duration)
            })

    # Create Table from list of dicts
    if not rows:
        return Table()

    # print progres
    print(f'Making table of {len(rows)} rows')

    t = Table(rows)

    # Convert the timestamp column to actual Astropy Time objects
    t['timestamp'] = Time(t['timestamp'], format='iso', scale='utc')
    t['runtime'].unit = 's'

    # Cumulative time based only on timestamp deltas between consecutive rows.
    cumulative_time = [0.0]
    for it in range(1, len(t)):
        dt_sec = (t['timestamp'][it] - t['timestamp'][it - 1]).to_value('sec')
        cumulative_time.append(cumulative_time[-1] + dt_sec)

    t['cumulative_time'] = cumulative_time
    t['cumulative_time'].unit = 's'

    return t

def plot_log(log_table, title=None, mask=None):


    # apply mask if given
    if mask is not None:
        plot_dates = log_table['timestamp'].plot_date

        mintime = Time(mask[0], format='iso').plot_date
        maxtime = Time(mask[1], format='iso').plot_date

        mask = (plot_dates >= mintime) & (plot_dates <= maxtime)
        log_table = log_table[mask]

    # print progress
    if title is not None:
        print(f'Plotting graph {title}')
    else:
        print('Plotting graph')

    timestamps_dt = log_table['timestamp'].to_datetime()
    # 3x2 layout: full-series plots on the left, 1-hour zoom panels on the right.
    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(nrows=3, ncols=2, width_ratios=[3, 1],
                          wspace=0.01, hspace=0.08)
    frames = [fig.add_subplot(gs[row, 0]) for row in range(3)]
    zoom_frames = [fig.add_subplot(gs[row, 1]) for row in range(3)]

    plot_dates = log_table['timestamp'].plot_date
    t_mid = 0.5 * (plot_dates.min() + plot_dates.max())
    half_window_days = (ZOOM_IN_SIZE / 2 ) / 24.0
    zoom_start = t_mid - half_window_days
    zoom_end = t_mid + half_window_days
    zoom_mask = (plot_dates >= zoom_start) & (plot_dates <= zoom_end)

    # work out the average time
    average = np.nanmean(log_table['runtime'])
    total = np.nanmax(log_table['cumulative_time'] / 3600)
    # plot time
    frames[0].plot(timestamps_dt, log_table['runtime'],
                   color='purple', ls='None', marker='.',
                   label=f'Run time [average={average} s]', alpha=0.25)
    zoom_frames[0].plot(timestamps_dt[zoom_mask],
                        log_table['runtime'][zoom_mask],
                        color='purple', ls='None', marker='.',
                        alpha=0.5)
    # plot the cumulative time taken (from timestamp only)
    frames[1].plot(timestamps_dt, log_table['cumulative_time'] / 3600,
                   color='orange', ls='None', marker='.',
                   label=f'Cumulative Time [total={total} hr]',
                   alpha=0.25)
    zoom_frames[1].plot(timestamps_dt[zoom_mask],
                        (log_table['cumulative_time'][zoom_mask] / 3600),
                        color='orange', ls='None', marker='.',
                        alpha=0.5)
    # plot the dt between cumulative time taken
    frames[2].plot(timestamps_dt[:-1],
                   np.diff(log_table['cumulative_time']),
                   color='blue', ls='-', marker='.',
                   label='Cumulative Time', alpha=0.5)
    zoom_frames[2].plot(timestamps_dt[:-1][zoom_mask[:-1]],
                        np.diff(log_table['cumulative_time'])[zoom_mask[:-1]],
                        color='blue', ls='-', marker='.',
                        alpha=0.7)

    for frame in frames + zoom_frames:
        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        frame.xaxis.set_major_locator(locator)
        frame.xaxis.set_major_formatter(formatter)

    for frame in zoom_frames:
        frame.set_xlim(zoom_start, zoom_end)
        frame.yaxis.set_label_position('right')
        frame.yaxis.tick_right()
        frame.yaxis.set_ticks_position('right')
        frame.tick_params(axis='y', left=False, labelleft=False,
                          right=True, labelright=True)

    # axis labels
    frames[0].set_ylabel('Run Time (seconds)')
    frames[1].set_xlabel('Timestamp')
    frames[1].set_ylabel('Cumulative Time [hr]')
    frames[2].set_ylabel('Delta Time [s] (between recipe runs)')
    zoom_frames[0].set_ylabel('Run Time (seconds)', labelpad=4)
    zoom_frames[1].set_ylabel('Cumulative Time [hr]', labelpad=4)
    zoom_frames[2].set_ylabel('Delta Time [s]', labelpad=4)
    zoom_frames[1].set_xlabel('Zoomed Timestamp')

    zoom_frames[0].set_title(f'{ZOOM_IN_SIZE} hr zoom around midpoint')

    for frame in frames:
        frame.legend(loc=0)

    if title is not None:
        title = str(title)
    elif RECIPE is None:
        title = 'Recipe Runtimes and Cumulative Times'
    else:
        title = f'{RECIPE} Runtimes and Cumulative Time'

    plt.suptitle(title)
    plt.subplots_adjust(bottom=0.05, top=0.95, right=0.94, left=0.04)

    plt.show(block=True)
    plt.close()


# =============================================================================
# Start of main code
# =============================================================================
if __name__ == '__main__':

    for case in [1,2,3]:
        if case == 1:
            FILE_PATH = '/home/cook//alliance_timing/neil_alliance_pp_time.txt'
            # define the start date of the log (the messages are in HH:MM:SS.SS
            START_DATE = '2026-03-17'
            # define title
            TITLE = f'Neil Alliance pp times [{START_DATE}]'
            # cut down the data (to avoid big jumps)
            MASK = None
        elif case == 2:
            FILE_PATH = '/home/cook//alliance_timing/lison_alliance_pp_time.txt'
            # define the start date of the log (the messages are in HH:MM:SS.SS
            START_DATE = '2025-11-14'
            # define title
            TITLE = f'Lison Alliance pp times [{START_DATE}]'
            # cut down the data (to avoid big jumps)
            MASK = None
        elif case == 3:
            FILE_PATH = '/home/cook//alliance_timing/sctitan_udem_pp_time.txt'
            # define the start date of the log (the messages are in HH:MM:SS.SS
            START_DATE = '2024-12-06'
            # define title
            TITLE = f'SCTITAN  pp times [{START_DATE}]'
            # cut down the data (to avoid big jumps)
            MASK = ['2024-12-10 00:00:00', '2024-12-11 00:00:00']
        else:
            raise ValueError('Invalid case number')


        # get times
        log_table = extract_log_to_table(FILE_PATH, start_date=START_DATE,
                                         recipe=RECIPE)
        # plot log table
        plot_log(log_table, title=TITLE, mask=MASK)

# =============================================================================
# End of main code
# =============================================================================
