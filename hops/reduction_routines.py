

import time
import numpy as np
import shutil
import hops.pylightcurve3 as plc
import matplotlib.cm as cm
import matplotlib.patches as mpatches

from astropy.io import fits as pf

from hops.hops_tools.windows import *
from hops.hops_tools.logs import log
from hops.hops_tools.fits import *

import glob


class Reduction:

    def __init__(self, progress_window):
        print('Reduction...')

        # get variables

        self.observation_files = find_fits_files(log.read_local_log('pipeline', 'observation_files'))
        self.bias_files = find_fits_files(log.read_local_log('reduction', 'bias_files'))
        self.dark_files = find_fits_files(log.read_local_log('reduction', 'dark_files'))
        self.flat_files = find_fits_files(log.read_local_log('reduction', 'flat_files'))

        self.reduction_directory = log.read_local_log('pipeline', 'reduction_directory')
        self.reduction_prefix = log.read_local_log('pipeline', 'reduction_prefix')
        self.exposure_time_key = log.read_local_log('pipeline_keywords', 'exposure_time_key')
        self.mean_key = log.read_local_log('pipeline_keywords', 'mean_key')
        self.std_key = log.read_local_log('pipeline_keywords', 'std_key')
        self.observation_date_key = log.read_local_log('pipeline_keywords', 'observation_date_key')
        self.observation_time_key = log.read_local_log('pipeline_keywords', 'observation_time_key')
        self.frame_low_std = log.read_local_log('windows', 'frame_low_std')
        self.frame_upper_std = log.read_local_log('windows', 'frame_upper_std')

        self.bin_fits = int(log.read_local_log('reduction', 'bin_fits'))
        self.bin_to = int(log.read_local_log('reduction', 'bin_to'))
        self.master_bias_method = log.read_local_log('reduction', 'master_bias_method')
        self.master_dark_method = log.read_local_log('reduction', 'master_dark_method')
        self.master_flat_method = log.read_local_log('reduction', 'master_flat_method')

        # check if reduction directory exists

        if os.path.isdir(self.reduction_directory):
            shutil.rmtree(self.reduction_directory)

        os.mkdir(self.reduction_directory)

        self.progress_window = progress_window

    def run(self):

        fits = get_fits_data(self.observation_files[0])

        mean = np.median(fits[0].data)
        std = plc.mad(fits[0].data) * 1.5
        self.progress_window.ax.imshow(fits[0].data, origin='lower', cmap=cm.Greys_r,
                  vmin=mean + self.frame_low_std * std,
                  vmax=mean + self.frame_upper_std * std)

        self.progress_window.canvas.draw()

        self.progress_window.show()

        master_bias = self.get_master_bias()
        if not self.progress_window.exit:
            master_dark = self.get_master_dark(master_bias)
            if not self.progress_window.exit:
                master_flat = self.get_master_flat(master_bias, master_dark)
                if not self.progress_window.exit:
                    x, y, z = self.reduce(master_bias, master_dark, master_flat)
                    if not self.progress_window.exit:
                        self.show_sky(x, y, z)
                        if not self.progress_window.exit:
                            log.write_local_log('pipeline', True, 'reduction_complete')

        self.progress_window.hide()
        self.progress_window.exit = False

    def get_master_bias(self):

        bias_frames = []
        percent = 0
        lt0 = time.time()
        for counter, bias_file in enumerate(self.bias_files):

            if self.progress_window.exit:
                return None

            fits = get_fits_data(bias_file)

            bias_frames.append(np.ones_like(fits[0].data) * fits[0].data)

            new_percent = round(100 * (counter + 1) / float(len(self.bias_files)), 1)
            if new_percent != percent:
                lt1 = time.time()
                rm_time = (100 - new_percent) * (lt1 - lt0) / new_percent
                hours = rm_time / 3600.0
                minutes = (hours - int(hours)) * 60
                seconds = (minutes - int(minutes)) * 60

                self.progress_window.progress_bar_1['value'] = new_percent
                self.progress_window.percent_label_1.configure(text='{0} % ({1}h {2}m {3}s left)'.format(new_percent, int(hours),
                                                                                    int(minutes), int(seconds)))

                percent = new_percent

                self.progress_window.progress_bar_1.update()
                self.progress_window.percent_label_1.update()

        if len(bias_frames) > 0:
            if self.master_bias_method == 'median':
                master_bias = np.median(bias_frames, 0)
            elif self.master_bias_method == 'mean':
                master_bias = np.mean(bias_frames, 0)
            else:
                master_bias = np.median(bias_frames, 0)
        else:
            master_bias = 0.0
            self.progress_window.progress_bar_1['value'] = 0
            self.progress_window.percent_label_1.configure(text='No bias frames')
            self.progress_window.progress_bar_1.update()
            self.progress_window.percent_label_1.update()

        print('Median Bias: ', round(np.median(master_bias), 3))

        return master_bias

    def get_master_dark(self, master_bias):

        dark_frames = []
        percent = 0
        lt0 = time.time()

        for counter, dark_file in enumerate(self.dark_files):

            if self.progress_window.exit:
                return None

            fits = get_fits_data(dark_file)

            dark_frame = np.ones_like(fits[0].data) * fits[0].data
            dark_frames.append((dark_frame - master_bias) / fits[0].header[self.exposure_time_key])

            new_percent = round(100 * (counter + 1) / float(len(self.dark_files)), 1)
            if new_percent != percent:
                lt1 = time.time()
                rm_time = (100 - new_percent) * (lt1 - lt0) / new_percent
                hours = rm_time / 3600.0
                minutes = (hours - int(hours)) * 60
                seconds = (minutes - int(minutes)) * 60

                self.progress_window.progress_bar_2['value'] = new_percent
                self.progress_window.percent_label_2.configure(text='{0} % ({1}h {2}m {3}s left)'.format(new_percent, int(hours),
                                                                                    int(minutes), int(seconds)))

                percent = new_percent

                self.progress_window.progress_bar_2.update()
                self.progress_window.percent_label_2.update()

        if len(dark_frames) > 0:
            if self.master_dark_method == 'median':
                master_dark = np.median(dark_frames, 0)
            elif self.master_dark_method == 'mean':
                master_dark = np.mean(dark_frames, 0)
            else:
                master_dark = np.median(dark_frames, 0)
        else:
            master_dark = 0.0
            self.progress_window.progress_bar_2['value'] = 0
            self.progress_window.percent_label_2.configure(text='No dark frames')
            self.progress_window.progress_bar_2.update()
            self.progress_window.percent_label_2.update()

        print('Median Dark: ', round(np.median(master_dark), 3))

        return master_dark

    def get_master_flat(self, master_bias, master_dark):

        flat_frames = []
        percent = 0
        lt0 = time.time()
        for counter, flat_file in enumerate(self.flat_files):

            if self.progress_window.exit:
                return None

            fits = get_fits_data(flat_file)

            flat_frame = np.ones_like(fits[0].data) * fits[0].data
            flat_frames.append(flat_frame - master_bias - fits[0].header[self.exposure_time_key] * master_dark)

            new_percent = round(100 * (counter + 1) / float(len(self.flat_files)), 1)
            if new_percent != percent:
                lt1 = time.time()
                rm_time = (100 - new_percent) * (lt1 - lt0) / new_percent
                hours = rm_time / 3600.0
                minutes = (hours - int(hours)) * 60
                seconds = (minutes - int(minutes)) * 60

                self.progress_window.progress_bar_3['value'] = new_percent
                self.progress_window.percent_label_3.configure(text='{0} % ({1}h {2}m {3}s left)'.format(new_percent, int(hours),
                                                                                    int(minutes), int(seconds)))

                percent = new_percent

                self.progress_window.progress_bar_3.update()
                self.progress_window.percent_label_3.update()

        print('Median of each Flat: ', ' '.join([str(round(np.median(ff))) for ff in flat_frames]))
        if len(flat_frames) > 0:
            if self.master_flat_method == 'mean':
                flat_frames = [ff / np.mean(ff) for ff in flat_frames]
                master_flat = np.mean(flat_frames, 0)
            else:
                flat_frames = [ff / np.median(ff) for ff in flat_frames]
                master_flat = np.median(flat_frames, 0)
            print('Median Flat: ', round(np.median(master_flat), 3))
            master_flat = master_flat / np.median(master_flat)
            master_flat = np.where(master_flat == 0, 1, master_flat)
        else:
            master_flat = 1.0
            self.progress_window.progress_bar_3['value'] = 0
            self.progress_window.percent_label_3.configure(text='No flat frames')
            self.progress_window.progress_bar_3.update()
            self.progress_window.percent_label_3.update()

        return master_flat

    def reduce(self, master_bias, master_dark, master_flat):

        # correct each observation_files file

        percent = 0
        lt0 = time.time()

        testx = []
        testy = []
        testz = []

        for counter, science_file in enumerate(self.observation_files):

            if self.progress_window.exit:
                return None, None, None

            self.progress_window.label_4.configure(text='Reducing data and calculating statistics: {0}'.format(os.path.split(science_file)[1]))
            self.progress_window.label_4.update()

            # correct it with master bias_files, master dark_files and master flat_files

            fits = get_fits_data(science_file)

            data_frame = np.ones_like(fits[0].data) * fits[0].data
            data_frame = (data_frame - master_bias - fits[0].header[self.exposure_time_key] * master_dark) / master_flat

            if self.bin_fits > 1:
                data_frame = plc.bin_frame(data_frame, self.bin_fits)

            try:
                distribution = plc.one_d_distribution(data_frame.flatten()[::int(200000.0/self.bin_to)], gaussian_fit=True)
                mean = distribution[2]
                std = distribution[3]
            except:
                mean = np.median(data_frame)
                std = plc.mad(data_frame) * 1.5

            if self.observation_date_key == self.observation_time_key:
                observation_time = ' '.join(fits[0].header[self.observation_date_key].split('T'))
            else:
                observation_time = ' '.join([fits[0].header[self.observation_date_key].split('T')[0],
                                             fits[0].header[self.observation_time_key]])

            observation_time = plc.UTC(observation_time)

            testx.append(observation_time.jd)
            testy.append(mean / fits[0].header[self.exposure_time_key])
            testz.append(std)

            fits[0].header.set(self.mean_key, mean)
            fits[0].header.set(self.std_key, std)

            # write the new fits file
            # important to keep it like this for windows!
            time_in_file = observation_time.utc.isoformat()
            time_in_file = time_in_file.split('.')[0]
            time_in_file = time_in_file.replace('-', '_').replace(':', '_').replace('T', '_')

            hdu = pf.ImageHDU(header=fits[0].header, data=np.array(data_frame, dtype=np.float32))

            plc.save_fits(pf.HDUList([pf.PrimaryHDU(), hdu]), '{0}{1}{2}{3}_{4}'.format(self.reduction_directory,
                                                 os.sep, self.reduction_prefix, time_in_file, science_file.split(os.sep)[-1]))

            if counter == 0:
                self.progress_window.ax.cla()
                self.progress_window.ax.imshow(data_frame[::2, ::2], origin='lower', cmap=cm.Greys_r,
                          vmin=fits[0].header[self.mean_key] + self.frame_low_std * fits[0].header[self.std_key],
                          vmax=fits[0].header[self.mean_key] + self.frame_upper_std * fits[0].header[self.std_key])
                self.progress_window.ax.axis('off')

                self.progress_window.canvas.draw()

            # counter

            new_percent = round(100 * (counter + 1) / float(len(self.observation_files)), 1)
            if new_percent != percent:
                lt1 = time.time()
                rm_time = (100 - new_percent) * (lt1 - lt0) / new_percent
                hours = rm_time / 3600.0
                minutes = (hours - int(hours)) * 60
                seconds = (minutes - int(minutes)) * 60

                self.progress_window.progress_bar_4['value'] = new_percent
                self.progress_window.percent_label_4.configure(text='{0} % ({1}h {2}m {3}s left)'.format(new_percent, int(hours),
                                                                                    int(minutes), int(seconds)))

                percent = new_percent

            self.progress_window.update()

        return testx, testy, testz
    #
    def show_sky(testx, testy, testz):
        pass
    #
    #     root = ProgressWindow('HOPS - Alignment', 0, 0, 5)
    #
    #     testx = np.array(np.array(testx) - testx[0]) * 24.0 * 60.0
    #
    #     reduction_trash_directory = log.read_local_log('pipeline', 'reduction_trash_directory')
    #     trash = log.read_local_log('pipeline', 'trash')
    #     if not trash:
    #         list_to_remove = []
    #     else:
    #         list_to_remove = list(np.int_(trash))
    #
    #     f = Figure()
    #     ax = f.add_subplot(2, 1, 1)
    #     ax2 = f.add_subplot(2, 2, 3)
    #     ax3 = f.add_subplot(2, 2, 4)
    #
    #     f.patch.set_facecolor('white')
    #     canvas = root.FigureCanvasTkAgg(f)
    #     canvas.get_tk_widget().pack()
    #     root.NavigationToolbar2Tk(canvas)
    #
    #     ax.plot(testx, testy, 'ko', ms=3)
    #     for ii in list_to_remove:
    #         ax.plot(testx[ii], testy[ii], 'ro', ms=3)
    #
    #     time_dt = np.median(np.array(testx[1:]) - np.array(testx[:-1]))
    #     arrow = mpatches.Arrow(testx[0], 0, 0, testy[0], width=time_dt, fc='r')
    #     ax.add_patch(arrow)
    #     ax.set_xlabel('Time (minutes from observation start)')
    #     ax.set_ylabel('Sky (counts / pixel / second)')
    #     ax.tick_params(top='on', bottom='off', labeltop='on', labelbottom='off')
    #     ax.xaxis.set_label_position('top')
    #
    #     science = find_fits_files(os.path.join(reduction_directory, '*'))
    #
    #     fits = pf.open(science[0], memmap=False)
    #
    #     ax2.imshow(fits[1].data, origin='lower', cmap=cm.Greys_r,
    #                vmin=fits[1].header[mean_key] + frame_low_std * fits[1].header[std_key],
    #                vmax=fits[1].header[mean_key] + frame_upper_std * fits[1].header[std_key])
    #     ax2.axis('off')
    #
    #     fits.close()
    #
    #     ax3.text(-100105, -100100, 'Select faulty frames', va='center', ha='center')
    #     ax3.text(-100111, -100101, '>On the time-sky graph above\n'
    #                                'double-click on a point to see\n'
    #                                'the frame on the left panel.\n'
    #                                '>To mark this point as faulty,\n'
    #                                'use the right double-click.\n'
    #                                '>To undo, use the right\n'
    #                                'double-click again.', va='top')
    #     ax3.text(-100105, -100109, 'RUN ALIGNMENT', bbox={'facecolor': 'white', 'alpha': 0.5, 'pad': 5},
    #              va='center', ha='center')
    #     ax3.set_xlim(-100110, -100100)
    #     ax3.set_ylim(-100110, -100100)
    #     ax3.axis('off')
    #
    #     def update_window_show(event):
    #
    #         if event.inaxes is not None:
    #
    #             if testx[0] < event.xdata < testx[-1]:
    #
    #                 plot_hjd = np.argmin(np.abs(event.xdata - np.array(testx)))
    #
    #                 if event.dblclick:
    #
    #                     del ax.patches[0]
    #                     arrow2 = mpatches.Arrow(testx[plot_hjd], 0, 0, testy[plot_hjd], width=time_dt, fc='r')
    #                     ax.add_patch(arrow2)
    #
    #                     ax2.cla()
    #                     fits2 = pf.open(science[plot_hjd], memmap=False)
    #                     ax2.imshow(fits2[1].data, origin='lower', cmap=cm.Greys_r,
    #                                vmin=fits2[1].header[mean_key] + frame_low_std * fits2[1].header[std_key],
    #                                vmax=fits2[1].header[mean_key] + frame_upper_std * fits2[1].header[std_key])
    #                     ax2.axis('off')
    #
    #                     if event.button == 3:
    #
    #                         pltxlim1, pltxlim2 = ax.get_xlim()
    #                         pltylim1, pltylim2 = ax.get_ylim()
    #
    #                         if plot_hjd not in list_to_remove:
    #                             ax.plot(testx[plot_hjd], testy[plot_hjd], 'ro', ms=3)
    #                             list_to_remove.append(plot_hjd)
    #                         else:
    #                             ax.plot(testx[plot_hjd], testy[plot_hjd], 'ko', ms=3)
    #                             list_to_remove.remove(plot_hjd)
    #
    #                         ax.set_xlim(pltxlim1, pltxlim2)
    #                         ax.set_ylim(pltylim1, pltylim2)
    #
    #                 canvas.draw()
    #
    #             else:
    #                 if -100110 < event.ydata < -100108:
    #                     if -100108 < event.xdata < -100102:
    #                         ax3.cla()
    #                         ax3.text(-100105, -100100, 'Select faulty frames', va='center', ha='center')
    #                         ax3.text(-100111, -100101, '>On the time-sky graph above\n'
    #                                                    'double-click on a point to see\n'
    #                                                    'the frame on the left panel.\n'
    #                                                    '>To mark this point as faulty,\n'
    #                                                    'use the right double-click.\n'
    #                                                    '>To undo, use the right\n'
    #                                                    'double-click again.', va='top')
    #                         ax3.text(-100105, -100109, 'RUN ALIGNMENT', color='w',
    #                                  bbox={'facecolor': 'blue', 'alpha': 0.5, 'pad': 5},
    #                                  va='center', ha='center')
    #                         ax3.set_xlim(-100110, -100100)
    #                         ax3.set_ylim(-100110, -100100)
    #                         ax3.axis('off')
    #                         canvas.draw()
    #                         time.sleep(0.5)
    #                         root.exit = True
    #
    #     f.canvas.callbacks.connect('button_press_event', update_window_show)
    #     # f.canvas.callbacks.connect('button_press_event', run_alignment)
    #
    #     root.show()
    #
    #     while not root.exit:
    #         root.update()
    #
    #     root.close()
    #     if not os.path.isdir(os.path.join(reduction_directory, reduction_trash_directory)):
    #         os.mkdir(os.path.join(reduction_directory, reduction_trash_directory))
    #
    #     for iii in list_to_remove:
    #         shutil.move(science[iii], os.path.join(reduction_directory, reduction_trash_directory))
    #
    #     log.write_local_log('pipeline', list(map(str, list_to_remove)), 'trash')
