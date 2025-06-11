#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-05-15 at 09:06

@author: cook
"""
import os
from sqlalchemy import create_engine, text
import time
import csv
from datetime import datetime
import argparse
import pandas as pd


# =============================================================================
# Define variables
# =============================================================================
# Configure your connection string
ENGINE_URL = "mysql+pymysql://{user}:{password}@{host}/{database}"

# List of InnoDB status keys to track
DB_KEYS = {"Innodb_buffer_pool_bytes_data",
           "Innodb_buffer_pool_pages_data",
           "Innodb_buffer_pool_pages_free",
           "Innodb_buffer_pool_read_requests",
           "Innodb_buffer_pool_reads",
           "Innodb_buffer_pool_write_requests",
           "Innodb_buffer_pool_pages_flushed"}

DB_SETUPS = dict()
DB_SETUPS['spirou'] = dict(user='spirou', password='Covid19!',
                           host='cosmos.astro.umontreal.ca', database='spirou')
DB_SETUPS['nirps'] = dict(user='nirps', password='Covid19!',
                          host='rali.astro.umontreal.ca', database='nirps')
# CSV file to write to
OUTPUT_FILE = "innodb_monitor_log.csv"

# =============================================================================
# Define functions
# =============================================================================
def write_csv(log_file, result):
    if os.path.exists(log_file):
        df = pd.read_csv(log_file)
    else:
        data = {key: [] for key in DB_KEYS}
        data["timestamp"] = []
        df = pd.DataFrame(data)

    # Filter for the keys we care about
    data = {row[0]: [int(row[1])] for row in result if row[0] in DB_KEYS}
    data["timestamp"] = [datetime.now().isoformat()]
    df_new = pd.DataFrame(data)

    df_save = pd.concat([df, df_new], ignore_index=True)
    df_save.to_csv(log_file, index=False)

    return df_new


def monitor(log_file, instrument):
    # get database arguments for this instrument
    db_args = DB_SETUPS[instrument]
    # create engine
    engine = create_engine(ENGINE_URL.format(**db_args))
    # start a counter
    counter = 0

    start = time.time()

    with engine.connect() as conn:

        while True:

            time_running = time.time() - start
            print('\n\n' + '*' * 50)
            print(f'* {counter + 1}  (Time running = {time_running:.4f})')
            print('*' * 50)

            result = conn.execute(text("SHOW GLOBAL STATUS")).fetchall()

            df_new = write_csv(log_file, result)

            print(df_new.iloc[0])
            time.sleep(5)  # Adjust the interval as needed

            # add to the counter
            counter += 1


class Monitor:
    def __init__(self, log_file):
        self.log_file = log_file
        self.fig = None
        self.frames = []
        self.lines = []
        self.xvalues = []
        self.yvalues = []
        self.data = dict()

        self.suptitle = ''


    def get_data(self):
        # Load CSV
        while True:
            try:
                df = pd.read_csv(self.log_file, parse_dates=["timestamp"])
                df.set_index("timestamp", inplace=True)
                break
            except Exception as _:
                print(f'Cannot load {self.log_file}. Trying again...')
                time.sleep(1)
        # Compute time difference between rows
        time_deltas = df.index.to_series().diff().dt.total_seconds()

        # Calculate raw diffs
        delta_df = df.diff().fillna(0)

        # Zero out rows where the time delta is too large (e.g., > 10s)
        invalid = time_deltas > 10  # or choose another threshold
        delta_df[invalid] = 0

        # Normalize by actual time difference (if not too large)
        time_deltas = time_deltas.where(time_deltas <= 10, 5)

        delta_df = delta_df.divide(time_deltas, axis=0).fillna(0)

        # Calculate hit ratio (avoid division by zero)
        ratio = df["Innodb_buffer_pool_reads"] / df["Innodb_buffer_pool_read_requests"]
        infs = [float('inf'), -float('inf')]
        df["hit_ratio"] = 1 - (ratio).replace(infs, 0).fillna(0)
        # push into storage
        self.data['x'] = df.index
        self.data['dx'] = delta_df.index
        self.data['l1a_y'] = delta_df["Innodb_buffer_pool_reads"]
        self.data['l1b_y'] = delta_df["Innodb_buffer_pool_read_requests"]
        self.data['l1c_y'] = delta_df["Innodb_buffer_pool_write_requests"]
        self.data['l1d_y'] = delta_df["Innodb_buffer_pool_pages_flushed"]
        self.data['l2_y'] = df["hit_ratio"]
        self.data['l3a_y'] = df["Innodb_buffer_pool_pages_data"]
        self.data['l3b_y'] = df["Innodb_buffer_pool_pages_free"]



    def plot(self):
        import matplotlib.pyplot as plt

        self.get_data()
        # plot intial plot
        plt.close()
        self.fig, frames = plt.subplots(nrows=3, ncols=1, figsize=(12, 6))
        self.frames = frames
        # ----------------------------
        # Plot 1: read/writes
        # ----------------------------
        line1a, = frames[0].plot(self.data['dx'], self.data['l1a_y'],
                                   label="Physical Reads/sec")
        line1b, = frames[0].plot(self.data['dx'], self.data['l1b_y'],
                                   label="Logical Reads/sec")
        line1c, = frames[0].plot(self.data['dx'], self.data['l1c_y'],
                                   label="Write Requests/sec")
        line1d, = frames[0].plot(self.data['dx'], self.data['l1d_y'],
                                   label="Pages Flushed/sec")
        self.lines += [line1a, line1b, line1c, line1d]
        self.xvalues += ['dx', 'dx', 'dx', 'dx']
        self.yvalues += ['l1a_y', 'l1b_y', 'l1c_y', 'l1d_y']

        frames[0].legend()
        frames[0].set_title("InnoDB I/O Activity Over Time")
        frames[0].set_ylabel("Ops/sec")
        frames[0].set_xlabel("Time")
        frames[0].grid(True)

        # ----------------------------
        # Plot 2: Cache Hit Ratio Over Time
        # ----------------------------
        line2, = frames[1].plot(self.data['x'], self.data['l2_y'],
                                  color='green', label="Cache Hit Ratio")
        frames[1].set_title("InnoDB Buffer Pool Cache Hit Ratio Over Time")
        frames[1].set_ylabel("Hit Ratio")
        frames[1].set_ylim(0, 1.05)
        frames[1].set_xlabel("Time")
        frames[1].grid(True)
        frames[1].legend()

        self.lines += [line2]
        self.xvalues += ['x']
        self.yvalues += ['l2_y']
        # ----------------------------
        # Plot 3: Buffer Pool Memory Usage
        # ----------------------------
        line3a, = frames[2].plot(self.data['x'], self.data['l3a_y'],
                                   label="Pages with Data", color='blue')
        line3b, = frames[2].plot(self.data['x'], self.data['l3b_y'],
                                   label="Free Pages", color='orange')
        frames[2].set_title("InnoDB Buffer Pool Memory Usage Over Time")
        frames[2].set_ylabel("Pages")
        frames[2].set_xlabel("Time")
        frames[2].grid(True)
        frames[2].legend()


        self.lines += [line3a, line3b]
        self.xvalues += ['x', 'x']
        self.yvalues += ['l3a_y', 'l3b_y']

        self.suptitle = plt.suptitle("Waiting for data...")

        plt.tight_layout()

    def update(self, counter):

        print(f'Counter = {counter}')

        if self.fig is None:
            return
        # re-get data
        self.get_data()

        for it in range(len(self.lines)):

            line = self.lines[it]
            xvalue = self.xvalues[it]
            yvalue = self.yvalues[it]

            line.set_xdata(self.data[xvalue])
            line.set_ydata(self.data[yvalue])

        for frame in self.frames:
            frame.relim()             # recompute axes limits
            frame.autoscale_view()    # autoscale
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

        self.suptitle.set_text(f'Counter = {counter}')

        time.sleep(1)  # pause before next update (adjust as needed)

    def run(self):
        import matplotlib.pyplot as plt

        plt.ion()
        # initial plot
        self.plot()

        counter = 0
        # interactive loop
        while True:
            try:
                self.update(counter)
            except Exception as e:
                print(e)
                time.sleep(1)
            counter += 1


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    # add a positional argument
    parser.add_argument('instrument', nargs=1, type=str,
                        choices=['nirps', 'spirou'],
                        help='Instrument to use')
    # optional argument for mode
    parser.add_argument('--mode', type=str, choices=['monitor', 'plot'],
                        default='monitor',
                        help='Mode to use "monitor" for getting stats or '
                             '"plot" for plotting the current stats')
    parser.add_argument('--log', type=str, default=OUTPUT_FILE,
                        help='Log file to write to')
    # get args
    args = parser.parse_args()
    # first argument is the instrument
    instrument = args.instrument[0]
    # get the log file to use
    output_file = args.log
    # -------------------------------------------------------------------------
    if args.mode == 'monitor':
        monitor(output_file, instrument)
    elif args.mode == 'plot':
        plotter = Monitor(output_file)
        plotter.run()

# =============================================================================
# End of code
# =============================================================================
