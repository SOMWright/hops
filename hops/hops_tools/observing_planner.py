
import datetime
import numpy as np
from astroquery.simbad import Simbad

import hops.pylightcurve3 as plc

from hops.hops_tools.windows import *
from hops.hops_tools.logs import log
from hops.hops_tools.tests import *


class ObservingPlanner(AddOnWindow):
    
    def __init__(self):
    
        AddOnWindow.__init__(self, 'Observing Planner')
    
        # get variables from log and set as tk variables those to be modified
    
        self.target = self.StringVar(' ')
        self.target_search = self.StringVar('Orion nebula')
        self.target_ra = self.StringVar(' ')
        self.target_dec = self.StringVar(' ')
        now = datetime.datetime.now()
        self.obs_year_month = self.StringVar('{0} {1}'.format(now.year, str(now.month).zfill(2)))
        self.latitude = self.StringVar(log.read_local_log_profile('observatory_lat'))
        self.longitude = self.StringVar(log.read_local_log_profile('observatory_long'))
        self.horizon_s = self.StringVar(log.read_local_log_profile('observatory_horizon_s'))
        self.horizon_sw = self.StringVar(log.read_local_log_profile('observatory_horizon_sw'))
        self.horizon_w = self.StringVar(log.read_local_log_profile('observatory_horizon_w'))
        self.horizon_nw = self.StringVar(log.read_local_log_profile('observatory_horizon_nw'))
        self.horizon_n = self.StringVar(log.read_local_log_profile('observatory_horizon_n'))
        self.horizon_ne = self.StringVar(log.read_local_log_profile('observatory_horizon_ne'))
        self.horizon_e = self.StringVar(log.read_local_log_profile('observatory_horizon_e'))
        self.horizon_se = self.StringVar(log.read_local_log_profile('observatory_horizon_se'))
        self.observatory = self.StringVar(log.read_local_log_profile('observatory'))
        self.timezone = self.IntVar(log.read_local_log_profile('observatory_time_zone'))

        # create widgets

        self.figure, self.canvas, figure_frame = self.Figure(figsize=(6, 6.5), show_nav=True)
        self.ax1 = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.1, right=1-0.05, bottom=0.1, top=0.85)

        self.target_ra_entry = self.Entry(textvariable=self.target_ra, width=30)
        self.target_ra_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.target_dec_entry = self.Entry(textvariable=self.target_dec, width=30)
        self.target_dec_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.target_ra_dec_test = self.StringVar(' ')

        self.obs_year_month_entry = self.Entry(textvariable=self.obs_year_month, width=30)
        self.obs_year_month_entry.bind(sequence='<KeyRelease>', func=self.update_window)
        self.obs_year_month_test = self.StringVar(' ')
    
        self.plot_button = self.Button(text='RESET PLOT', command=self.plot)

        # setup window

        self.setup_window([
            [[figure_frame, 4, 1, 30]],
            [[self.Label(text='Observatory'), 2]],
            [[self.Entry(textvariable=self.observatory, width=30), 2]],
            [],
            [[self.Label(text='Latitude'), 1],
             [self.Label(text='Horizon altitude S (deg)'), 2],
             [self.Entry(textvariable=self.horizon_s, width=10), 3]],
            [[self.Entry(textvariable=self.latitude), 1],
             [self.Label(text='Horizon altitude SW (deg)'), 2],
             [self.Entry(textvariable=self.horizon_sw, width=10), 3]],
            [[self.Label(text='Horizon altitude W (deg)'), 2],
             [self.Entry(textvariable=self.horizon_w, width=10), 3]],
            [[self.Label(text='Longitude'), 1],
             [self.Label(text='Horizon altitude NW (deg)'), 2],
             [self.Entry(textvariable=self.horizon_nw, width=10), 3]],
            [[self.Entry(textvariable=self.longitude), 1],
             [self.Label(text='Horizon altitude N (deg)'), 2],
             [self.Entry(textvariable=self.horizon_n, width=10), 3]],
            [[self.Label(text='Horizon altitude NE (deg)'), 2],
             [self.Entry(textvariable=self.horizon_ne, width=10), 3]],
            [[self.Label(text='Time Zone'), 1],
             [self.Label(text='Horizon altitude E (deg)'), 2],
             [self.Entry(textvariable=self.horizon_e, width=10), 3]],
            [[self.Entry(textvariable=self.timezone), 1],
             [self.Label(text='Horizon altitude SE (deg)'), 2],
             [self.Entry(textvariable=self.horizon_se, width=10), 3]],
            [],
            [[self.Label(text='Object'), 1]],
            [[self.Entry(textvariable=self.target_search), 1],
             [self.Button(text='SEARCH RA/DEC', command=self.search_object), 2], [self.Label(textvar=self.target), 3]],
            [],
            [[self.Label(text='Manual target RA DEC\n(hh:mm:ss +/-dd:mm:ss)'), 1, 1, 2], [self.target_ra_entry, 2],
             [self.Label(textvar=self.target_ra_dec_test), 3, 1, 2]],
            [[self.target_dec_entry, 2]],
            [[self.Label(text='Observation year and month\n(yyyy mm)'), 1], [self.obs_year_month_entry, 2],
             [self.Label(textvar=self.obs_year_month_test), 3]],
            [],
            [[self.plot_button, 2]],
            [],
        ])

    # define the function that updates the window

    def first_call(self):
        self.search_object()
        self.show()

    def update_window(self, *event):

        if not event:
            pass

        check_ra_dec = test_coordinates2(self.target_ra_entry.get(), self.target_dec_entry.get())
        self.target_ra_dec_test.set(check_ra_dec[1])
        check_year_month = test_date(self.obs_year_month_entry.get())
        self.obs_year_month_test.set(check_year_month[1])

        if check_ra_dec[0] and check_year_month[0]:
            self.plot_button['state'] = self.NORMAL
        else:
            self.plot_button['state'] = self.DISABLED

    def search_object(self):

        try:
            result_table = Simbad.query_object(self.target_search.get(), wildcard=True)[0]

            if result_table:
                self.target.set(result_table['MAIN_ID'].decode("utf-8"))
                self.target_ra.set(str(result_table['RA']))
                self.target_dec.set(str(result_table['DEC']))
            else:
                self.target_ra.set('')
                self.target_dec.set('')
                self.target.set('No results')

        except:
            self.target_ra.set('')
            self.target_dec.set('')
            self.target.set('No results')

        self.plot()
        self.update_window()

    def plot(self):
        self.ax1.cla()

        try:
            self.avc_plot(self.latitude.get(), self.longitude.get(), self.timezone.get(),
                          [
                              [0, self.horizon_s.get()],
                              [45, self.horizon_sw.get()],
                              [90, self.horizon_w.get()],
                              [135, self.horizon_nw.get()],
                              [180, self.horizon_n.get()],
                              [225, self.horizon_ne.get()],
                              [270, self.horizon_e.get()],
                              [315, self.horizon_se.get()],

                         ],
                         self.target_ra.get(), self.target_dec.get(), self.obs_year_month.get(), self.ax1,
                         self.target.get(), self.observatory.get())

        except:
            pass

        self.canvas.draw()

    def target_azimuth_altitude(self, target_ra, target_dec, observatory_latitude, sidereal_time):

        observatory_latitude *= np.pi / 180
        target_dec *= np.pi / 180

        ha = sidereal_time - target_ra / 15
        if ha < 0:
            ha += 24

        altitude = np.arcsin(np.clip(np.sin(target_dec) * np.sin(observatory_latitude)
                                     + np.cos(target_dec) * np.cos(observatory_latitude) * np.cos(
            ha * 15 * np.pi / 180), -1, 1))

        azimuth = np.pi - np.arccos(np.clip((np.sin(target_dec) - np.sin(altitude) * np.sin(observatory_latitude)) /
                                            (np.cos(altitude) * np.cos(observatory_latitude)), -1, 1))

        if ha >= 12:
            azimuth = 2 * np.pi - azimuth

        return azimuth * 180 / np.pi, altitude * 180 / np.pi

    def get_target_events(self, target_ra, target_dec, observatory_latitude, horizon):

        # horizon

        horizon_list = []
        for horizon_line in horizon.split('\n'):
            if horizon_line != '':
                horizon_list.append(horizon_line.split())
        if len(horizon_list) == 1:
            def horizon(azimuth):
                return float(horizon_list[0][0])
        else:
            horizon_list.append([360.0, horizon_list[0][0]])
            horizon_data = np.swapaxes(np.array(horizon_list, dtype=float), 0, 1)
            horizon = plc.interp1d(horizon_data[0], horizon_data[1])

        # horizon

        sidereal_time_list = list(np.arange(0, 24, 0.1)) + [24]
        sidereal_time_altitude_list = []

        for sidereal_time in sidereal_time_list:
            azimuth, altitude = self.target_azimuth_altitude(target_ra, target_dec, observatory_latitude, sidereal_time)
            sidereal_time_altitude_list.append([sidereal_time, altitude - horizon(azimuth)])

        sidereal_time_altitude_list.sort()
        sidereal_time_altitude_list = np.swapaxes(sidereal_time_altitude_list, 0, 1)
        sidereal_time_altitude_function = plc.interp1d(np.array(sidereal_time_altitude_list[0]),
                                                       np.array(sidereal_time_altitude_list[1]))

        def target_horizon_diference(x, st):
            return sidereal_time_altitude_function(st)

        # find events

        events = []

        test_alt = sidereal_time_altitude_function(sidereal_time_list[0])
        for sidereal_time in sidereal_time_list:
            new_test_alt = sidereal_time_altitude_function(sidereal_time)
            if test_alt * new_test_alt < 0:
                popt, pcov = plc.curve_fit(target_horizon_diference, [0], [0], p0=[sidereal_time - 0.1])
                if test_alt < new_test_alt:
                    events.append([popt[0], 'target_rise'])
                else:
                    events.append([popt[0], 'target_set'])
            else:
                pass

            test_alt = new_test_alt

        return events, sidereal_time_altitude_function

    def avc_plot(self, latitude, longitude, tmzn, horizon, target_ra, target_dec, year_mont_string, ax, name,
                 observatory_name):
        ax.cla()

        target = plc.Target(plc.Hours(target_ra), plc.Degrees(target_dec))
        observatory = plc.Observatory(plc.Degrees(latitude), plc.Degrees(longitude), tmzn, horizon)
        observation = plc.Observation(target, observatory)

        year = int(year_mont_string.split()[0])
        month = int(year_mont_string.split()[1])

        months = ['xx', 'January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September',
                  'October',
                  'November', 'December']
        if (year - 2000) / 4.0 - int((year - 2000) / 4.0) == 0:
            days = [0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        else:
            days = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        time_0 = plc.LT('{0}-{1}-1 12:00:00'.format(year, month), observatory)
        jd_0 = time_0.jd

        time_1 = plc.JD(jd_0 + days[month])

        events = []

        # mid-day splits

        for jj in range(days[month] + 1):
            events.append([plc.JD(time_0.jd + jj), 'mid-day'])

        # target rise/set events

        events += observation.rise_set_events(time_0, time_1)

        # sun rise/set events

        for jj in range(days[month]):

            check_time = plc.JD(time_0.jd + jj)
            check_st = check_time.lst(observatory).hours

            sun = check_time.get_sun()

            # sun rise/set

            if - (90 - observatory.latitude.deg_pm) < sun.dec.deg_pm < 90 - observatory.latitude.deg_pm:

                rise_ha = np.arccos(- sun.dec.tan * observatory.latitude.tan) * 12 / np.pi

                if rise_ha < 12:
                    set_ha = rise_ha
                    rise_ha = 24 - rise_ha
                else:
                    set_ha = 24 - rise_ha

                rise_st = rise_ha + sun.ra.hours
                if rise_st > 24:
                    rise_st -= 24

                set_st = set_ha + sun.ra.hours
                if set_st > 24:
                    set_st -= 24

                if rise_st < check_st:
                    next_rise_in_st_hours = 24 + rise_st - check_st
                else:
                    next_rise_in_st_hours = rise_st - check_st

                if set_st < check_st:
                    next_set_in_st_hours = 24 + set_st - check_st
                else:
                    next_set_in_st_hours = set_st - check_st

                dt = next_rise_in_st_hours * (23.9344696 / 24)
                if dt < 24:
                    events.append([plc.JD(check_time.jd + dt / 24), 'sun_rise'])

                dt = next_set_in_st_hours * (23.9344696 / 24)
                if dt < 24:
                    events.append([plc.JD(check_time.jd + dt / 24), 'sun_set'])

            # sun -18 rise/set

            if - (90 - observatory.latitude.deg_pm + 18.0) < sun.dec.deg_pm < 90 - (observatory.latitude.deg_pm + 18):

                rise_ha = np.arccos(np.sin((-18.0) * np.pi / 180) / sun.dec.cos / observatory.latitude.cos
                                    - sun.dec.tan * observatory.latitude.tan) * 12 / np.pi
                if rise_ha < 12:
                    set_ha = rise_ha
                    rise_ha = 24 - rise_ha
                else:
                    set_ha = 24 - rise_ha

                rise_st = rise_ha + sun.ra.hours
                if rise_st > 24:
                    rise_st -= 24

                set_st = set_ha + sun.ra.hours
                if set_st > 24:
                    set_st -= 24

                if rise_st < check_st:
                    next_rise_in_st_hours = 24 + rise_st - check_st
                else:
                    next_rise_in_st_hours = rise_st - check_st

                if set_st < check_st:
                    next_set_in_st_hours = 24 + set_st - check_st
                else:
                    next_set_in_st_hours = set_st - check_st

                dt = next_rise_in_st_hours * (23.9344696 / 24)
                if dt < 24:
                    events.append([plc.JD(check_time.jd + dt / 24), 'sun_rise_18'])

                dt = next_set_in_st_hours * (23.9344696 / 24)
                if dt < 24:
                    events.append([plc.JD(check_time.jd + dt / 24), 'sun_set_18'])

        events2 = [[ff[0].jd, ff[0], ff[1]] for ff in events]
        events2.sort(key=lambda ff: ff[0])

        #
        maxalt = str(round(observation.max_altitude.deg_pm, 1))

        ax.xaxis.tick_top()

        ax.set_title(observatory_name + '\n' + name + '   ' + months[month] + ' ' + str(year) +
                     '    max. alt. = ' + maxalt + ' degrees')
        ax.set_xlim((0, 1))
        ax.set_xlabel('HOUR (UTC{0:+.1f})'.format(tmzn))
        ax.set_xticks(np.arange(0, 24.5, 1))
        ax.set_xticklabels(('12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23',
                            '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'))
        ax.set_ylim((0, days[month] + 1))
        ax.set_ylabel('DAY')
        ax.set_yticks(np.arange(1, days[month] + 0.5, 1))
        ax.tick_params(bottom=True, top=True, left=True, right=True, labelbottom=True, labeltop=True,
                       labelright=False, labelleft=True)
        ax.grid(True, axis='y', linestyle='--')

        check_full_moon = plc.UTC('2000-1-21 04:41:00')

        for jj, ii in enumerate(events2[:-1]):

            moonphase = (np.sin((float(ii[0] + 0.5) - float(check_full_moon.jd)) * np.pi / 29.530589)) ** 2

            test_jd = 0.5 * (ii[0] + events2[jj + 1][0])
            dt_jd = 0.5 * (events2[jj + 1][0] - ii[0])

            day = 1 + int(test_jd - jd_0)

            time_range = [(ii[0] - jd_0 - int(ii[0] - jd_0)) * 24,
                          (events2[jj + 1][0] - jd_0 - int(events2[jj + 1][0] - jd_0)) * 24]
            if time_range[1] == 0:
                time_range[1] = 24

            alpha = 1
            if not observation.is_target_visible(plc.JD(ii[0] + dt_jd)):
                color = 'w'
                alpha = 0
            else:
                sun_az, sun_alt = observation.sun_azimuth_altitude(plc.JD(ii[0] + dt_jd))
                if sun_alt.deg_pm > 0:
                    color = 'y'
                elif sun_alt.deg_pm > -18:
                    color = 'r'
                else:
                    color = str(0.8 * (1 - moonphase))

            ax.plot(time_range, [day, day], linewidth=2.5, color=color, alpha=alpha)

            shift = {'left': +0.3, 'right': -0.3}

            if ii[2] == 'target_set':
                ax.plot(time_range[0], day, 'k*', mec='k', markersize=8)
                if time_range[0] > 20.5:
                    align = 'right'
                else:
                    align = 'left'
                ax.text(time_range[0] + shift[align], day + 0.4,
                        'set: ' + (ii[1].utc + datetime.timedelta(days=tmzn / 24)).isoformat().split('T')[1][:5],
                        va='center', ha=align, fontsize=9)

            if ii[2] == 'target_rise':
                ax.plot(time_range[0], day, 'w*', mec='k', markersize=8, markeredgewidth=0.5)
                if time_range[0] < 3.5:
                    align = 'left'
                else:
                    align = 'right'
                ax.text(time_range[0] + shift[align], day + 0.4,
                        'rise: ' + (ii[1].utc + datetime.timedelta(days=tmzn / 24)).isoformat().split('T')[1][:5],
                        va='center', ha=align, fontsize=9)
