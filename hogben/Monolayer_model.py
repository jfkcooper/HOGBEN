import os

import refnx.dataset
import refnx.reflect
import refnx.analysis

from hogben.models.base import BaseLipid

class MonolayerEllen(BaseLipid):
    """Defines a model describing a monolayer.

    Attributes:
        name (str): name of the monolayer sample.
        labels (list): label for each measured contrast.
        distances (numpy.ndarray): SLD profile x-axis range.
        deuterated (bool): whether the tailgroups are deuterated or not.
        scales (list): experimental scale factor for each measured contrast.
        bkgs (list): level of instrument background noise for each contrast.
        dq (float): instrument resolution.
        air_tg_rough (refnx.analysis.Parameter): air/tailgroup roughness.
        lipid_apm (refnx.analysis.Parameter): lipid area per molecule.
        hg_waters (refnx.analysis.Parameter): amount of headgroup bound water.
        monolayer_rough (refnx.analysis.Parameter): monolayer roughness.
        lipid_vf (refnx.analysis.Parameter): lipid hydration.
        params (list): varying model parameters.

    """

    def __init__(self, deuterated: bool = False):
        """Initialise a monolayer sample and set the sample properties"""
        self.name = 'Some_monolayer'


        self.labels = ['hDPPG-D2O', 'dDPPG-NRW', 'hDPPG-NRW']
        self.distances = None  # np.linspace(-25, 90, 500)
        self.deuterated = deuterated

        self.scale = [1.8899, 1.8832, 1.8574]
        self.bkg = [3.565e-6, 5.348e-6, 6.542e-6]
        self.dq = [3, 3, 3]

        # Define the varying parameters of the model.
        self.air_tg_rough = refnx.analysis.Parameter(
            5.0000, 'Air/Tailgroup Roughness', (5, 8)
        )
        self.lipid_apm = refnx.analysis.Parameter(
            54.1039, 'Lipid Area Per Molecule', (30, 80)
        )
        self.hg_waters = refnx.analysis.Parameter(
            6.6874, 'Headgroup Bound Waters', (0, 20)
        )
        self.monolayer_rough = refnx.analysis.Parameter(
            2.0233, 'Monolayer Roughness', (0, 10)
        )


        # We are only interested in the lipid APM.
        self.params = [self.lipid_apm]

        # Vary all of the parameters defined above.
        for param in self.params:
            param.vary = True

        # Call the BaseLipid constructor.
        super().__init__()

    def _create_objectives(self) -> None:
        """Creates objectives corresponding to each measured contrast.

        """
        # Define contrasts for below structures
        h2o, nrw, d2o = -0.5, 0.1, 6.35
        
        self.structures = [
        self._using_conditions(d2o),
        self._using_conditions(h2o),
        self._using_conditions(nrw),
        ]
        self.objectives = []



    def _using_conditions(self,
                          contrast_sld: float,
                          underlayers: Optional[list] = None,
                          deuterated: Optional[bool] = None,
                          ) -> refnx.reflect.Structure:
        """Creates a structure representing the monolayer measured using
           given measurement conditions.

        Args:
            contrast_sld (float): SLD of contrast to simulate.
            underlayers (list): thickness and SLD of each underlayer to add.
            deuterated (bool): whether the tailgroups are deuterated or not.
        Returns:
            refnx.reflect.Structure: monolayer using measurement conditions.

        """
        # Underlayers are not supported for this model.
        assert underlayers is None

        # If not specified, use the object's attribute.
        if deuterated is None:
            deuterated = self.deuterated

        # Convert the units of the given SLD.
        contrast_sld *= 1e-6

        # Define known SLDs of deuterated and hydrogenated
        d2o_sld = 6.35e-6
        h2o_sld = -0.56e-6

        # Define neutron scattering lengths.
        carbon_sl = 0.6646e-4
        oxygen_sl = 0.5843e-4
        hydrogen_sl = -0.3739e-4
        phosphorus_sl = 0.5130e-4
        deuterium_sl = 0.6671e-4

        # Calculate the total scattering length in each fragment.
        ch2_sl = 1 * carbon_sl + 2 * hydrogen_sl
        ch3_sl = 1 * carbon_sl + 3 * hydrogen_sl
        co2_sl = 1 * carbon_sl + 2 * oxygen_sl
        c3h5_sl = 3 * carbon_sl + 5 * hydrogen_sl
        po4_sl = 1 * phosphorus_sl + 4 * oxygen_sl
        h2o_sl = 2 * hydrogen_sl + 1 * oxygen_sl
        d2o_sl = 2 * deuterium_sl + 1 * oxygen_sl
        cd2_sl = 1 * carbon_sl + 2 * deuterium_sl
        cd3_sl = 1 * carbon_sl + 3 * deuterium_sl

        # Volumes of each fragment.
        ch2_vol = 28.1
        ch3_vol = 52.7 / 2
        co2_vol = 39.0
        c3h5_vol = 68.8
        po4_vol = 53.7
        water_vol = 30.4

        # Calculate volumes from components.
        hg_vol = po4_vol + 2 * c3h5_vol + 2 * co2_vol
        tg_vol = 28 * ch2_vol + 2 * ch3_vol

        # Calculate mole fraction of D2O from the bulk SLD.
        x = (contrast_sld - h2o_sld) / (d2o_sld - h2o_sld)

        # Calculate 'average' scattering length sum per water molecule.
        sl_sum_water = x * d2o_sl + (1 - x) * h2o_sl

        # Calculate scattering length sums for the other fragments.
        sl_sum_hg = po4_sl + 2 * c3h5_sl + 2 * co2_sl
        sl_sum_tg_d = 28 * cd2_sl + 2 * cd3_sl

        # Need to include number of hydrating water molecules in headgroup.
        hg_vol += water_vol * self.hg_waters
        sl_sum_hg += sl_sum_water * self.hg_waters

        hg_sld = sl_sum_hg / hg_vol

        # Calculate the SLD of the deuterated tailgroups.
        tg_d_sld = sl_sum_tg_d / tg_vol


        # Tailgroup and headgroup thicknesses.
        tg_thick = tg_vol / self.lipid_apm
        hg_thick = hg_vol / self.lipid_apm

        # Define the structure.
        air = refnx.reflect.SLD(0, name='Air')
        tg_d = refnx.reflect.SLD(tg_d_sld * 1e6, name='Deuterated Tailgroup')(
            tg_thick, self.air_tg_rough
        )
        hg = refnx.reflect.SLD(hg_sld * 1e6, name='Headgroup')(
            hg_thick, self.monolayer_rough
        )
        water = refnx.reflect.SLD(contrast_sld * 1e6, name='Water')(
            0, self.air_tg_rough
        )

        return air | tg_d | hg | water

    def sld_profile(self, save_path: str) -> None:
        """Plots the SLD profiles of the monolayer sample.

        Args:
            save_path (str): path to directory to save SLD profile to.

        """
        # Plot the SLD profile
        self._create_objectives()
        super().sld_profile(save_path, 'sld_profile',
                            ylim=(-0.6, 7.5))


