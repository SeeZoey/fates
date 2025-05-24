"""
Concrete class for running the patchprofile functional tests for FATES.
"""
import os
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from utils import round_up, get_color_palette, blank_plot
from functional_class import FunctionalTest



class PatchTest(FunctionalTest):
    """Patch test class
    """

    name = "patch"

    def __init__(self, test_dict):
        super().__init__(
            PatchTest.name,
            test_dict["test_dir"],
            test_dict["test_exe"],
            test_dict["out_file"],
            test_dict["use_param_file"],
            test_dict["other_args"],
        )
        self.plot = True

    def plot_output(self, run_dir: str, save_figs: bool, plot_dir: str):
        """Plots all patchprofile plots

        Args:
            run_dir (str): run directory
            out_file (str): output file name
            save_figs (bool): whether or not to save the figures
            plot_dir (str): plot directory to save the figures to
        """

        # read in patchprofile data
        patchprofile_dat = xr.open_dataset(os.path.join(run_dir, self.out_file))

        self.plot_lai_profile(patchprofile_dat, save_figs, plot_dir)

    @staticmethod
    def plot_lai_profile(data: xr.Dataset, save_fig: bool, plot_dir: str = None):
        """
        Plot average LAI profile over height, from a patchprofile NetCDF output

        Args:
            data (xr.Dataset): opened xarray dataset
            save_fig (bool): whether to save figure
            plot_dir (str): folder to save figure
        """

        import matplotlib.pyplot as plt
        import string
        import numpy as np
        import pandas as pd

        # read variables
        values = data.idspft.values
        ipft   = values[(values > 0) & (values < 10000)]
        dbh    = data.dbh.values[0:len(ipft)]
        tlai   = data.tlai_profile.values
        heightb = data.height_bottom.values

        # extract effective PFTs
        pft_indices = np.unique(ipft - 1)
        tlai_profile = tlai[:, pft_indices, :]
        heightb = heightb[:, pft_indices, :]

        flat_height = heightb.flatten()
        flat_lai = tlai_profile.flatten()

        # ranking the heights and lai
        df = pd.DataFrame({'height': flat_height, 'lai': flat_lai})
        df = df.groupby('height', as_index=False).mean()
        df = df.sort_values('height', ascending=False)

        plt.rcParams['font.family'] = 'Roboto'
        plt.rc('font', size=20)
        fig = plt.figure(figsize=(10, 12), dpi=300)
        ax1 = fig.add_axes([0.15, 0.1, 0.7, 0.8])  # [left, bottom, width, height]

        ax1.scatter(df["lai"], df["height"], c='r', marker='o', s=50, alpha=0.75)

        ax1.set_xlabel("LAI (m$^2$/m$^2$)", fontsize=20)
        ax1.set_ylabel("Bottom height of leaf layer (m)", fontsize=20)
        ax1.set_xlim(0, df["lai"].max()*1.1)

        ax1.tick_params(axis='both', which='major', width=2, length=8, labelsize=18)
        ax1.tick_params(axis='both', which='minor', width=1.5, length=5)
        [x.set_linewidth(2) for x in ax1.spines.values()]
        [x.set_color('k') for x in ax1.spines.values()]


        ax1.set_title(f'DBH = {dbh[:]}, PFT = {ipft[:]}', fontsize=18)

        if save_fig:
            os.makedirs(plot_dir, exist_ok=True)
            plt.savefig(os.path.join(plot_dir, "lai_profile.png"))
        else:
            plt.show()


