#!/usr/bin/env pythno
import numpy as np
import lusee
from pyuvdata.uvbeam import UVBeam

antenna_sim_path = "../../../AntennaSimResults/"
fname = "004_Freq1-50MHz_Delta1MHz_AntennaLength1-6m_Delta1m_AntennaAngle75deg_LanderHeight2m/RadiatedElectricField_AntennaLength6m_AntennaAngle75deg_LanderHeight2m_LBoxZ70cm_monopole_Phase+0deg.fits"

B = lusee.LBeam(antenna_sim_path+'/'+fname)
Nfreq, Ntheta, Nphi, Naxes = B.E.shape

uvb = UVBeam()
uvb.Naxes_vec = Naxes
uvb.Nfreqs = Nfreq
uvb.freq_array = np.ones((1,Nfreq))
uvb.freq_array[0,:] = B.freq*1e6
uvb.bandpass_array = np.ones((1,Nfreq))
uvb.Nspws = 1
uvb.spw_array = np.array([0])
uvb.antenna_type = "simple"
uvb.beam_type = 'efield'
uvb.pixel_coordinate_system = "az_za"
uvb.Naxes1 = Nphi
uvb.axis1_array = B.phi
uvb.Naxes2 = Ntheta
uvb.axis2_array = B.theta
uvb.Nfeeds = 1
data_shape = (uvb.Naxes_vec, uvb.Nspws, uvb.Nfeeds,
              uvb.Nfreqs, uvb.Naxes2, uvb.Naxes1)
uvb.data_array = np.zeros(data_shape,complex)
uvb.data_array [:,0,0,:,:,:] = B.E.transpose((3,0,1,2))
uvb.data_normalization = "physical" # unsure about this
uvb.telescope_name = "LuSEE-Night"
uvb.history = "test"
uvb.model_name = "Kaja 004"
uvb.model_version = "v0.1"
uvb.feed_version= "v0.1"
uvb.feed_array=["N"]
for ofs,c in enumerate([x for x in "NESW"]):
    uvb.feed_name = "LuSEE North  "+c
    cB = B.rotate(-90*ofs)
    uvb.axis1_array=cB.phi
    uvb.data_array [:,0,0,:,:,:] = cB.E.transpose((3,0,1,2))
    uvb.check()
    uvb.write_beamfits(f"data/lusee_004_{c}.fits", clobber=True)




