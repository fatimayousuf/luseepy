#!/usr/bin/env pythno
import numpy as np
import lusee
import pyuvdata
import pickle
from pyuvdata.uvbeam import UVBeam

print (pyuvdata.__version__)
#antenna_sim_path = "../../../../lusee_sky_simulations/AntennaSimResults"
#fname = "004_Freq1-50MHz_Delta1MHz_AntennaLength1-6m_Delta1m_AntennaAngle75deg_LanderHeight2m/RadiatedElectricField_AntennaLength6m_AntennaAngle75deg_LanderHeight2m_LBoxZ70cm_monopole_Phase+0deg.fits"

antenna_sim_path = "../../../AntennaSimResults/"
fname = "004_Freq1-50MHz_Delta1MHz_AntennaLength1-6m_Delta1m_AntennaAngle75deg_LanderHeight2m/RadiatedElectricField_AntennaLength6m_AntennaAngle75deg_LanderHeight2m_LBoxZ70cm_monopole_Phase+0deg.fits"


B = lusee.LBeam(antenna_sim_path+'/'+fname)
print (B.theta_deg)
Eproj = B.project_to_phi_theta()
Nfreq, Ntheta, Nphi, Naxes = Eproj.shape

uvb = UVBeam()
uvb.Naxes_vec = Naxes
uvb.Ncomponents_vec = 2
uvb.Nfreqs = Nfreq
uvb.freq_array = np.ones((1,Nfreq))
uvb.freq_array[0,:] = B.freq*1e6
uvb.bandpass_array = np.ones((1,Nfreq))
uvb.Nspws = 1
uvb.spw_array = np.array([0])
uvb.antenna_type = "simple"
#uvb.beam_type = 'efield'
uvb._set_efield()
uvb.pixel_coordinate_system = "az_za"
uvb._set_cs_params()
uvb.Naxes1 = Nphi
uvb.axis1_array = B.phi
uvb.Naxes2 = Ntheta
uvb.axis2_array = B.theta
uvb.Nfeeds = 2
uvb.feed_name = "LuSEE"
data_shape = (uvb.Naxes_vec, uvb.Nspws, uvb.Nfeeds,
              uvb.Nfreqs, uvb.Naxes2, uvb.Naxes1)
uvb.data_array = np.zeros(data_shape,complex)

uvb.data_normalization = "physical" # unsure about this
uvb.telescope_name = "LuSEE-Night"
uvb.history = "test"
uvb.model_name = "Kaja 004"
uvb.model_version = "v0.1"
uvb.feed_version= "v0.1"
uvb.feed_array=["X","Y"]
uvb.basis_vector_array = np.zeros(
        (uvb.Naxes_vec, uvb.Ncomponents_vec, uvb.Naxes2, uvb.Naxes1)
        )
uvb.basis_vector_array[0, 0, :, :] = 1.0
uvb.basis_vector_array[1, 1, :, :] = 1.0
for ofs,c in enumerate(["NE","SW"]):

    cB = B.rotate(-90*(2*ofs))
    cE= cB.project_to_phi_theta()
    uvb.data_array [:,0,0,:,:,:] = cE.transpose((3,0,1,2))
    cB = B.rotate(-90*(2*ofs+1))
    cE= cB.project_to_phi_theta()
    uvb.data_array [:,0,1,:,:,:] = cE.transpose((3,0,1,2))
    #uvb.data_array *= np.exp(3j)
    for i, theta in enumerate(B.theta):
        print(theta/np.pi*180, np.exp(-3*theta**2))
        uvb.data_array[:,0,0,:,i,:]*=np.exp(-4*theta**2)
        #uvb.data_array[:,0,1,:,i,:]*=np.exp(-3*theta**2)
    #uvb.data_array[:,:,:,:,:,(B.phi>np.pi)]=0
    

    #uvb.data_array*=10
    
    uvb.check()
    fname = f"data/lusee_004_{c}.fits"
    uvb.write_beamfits(fname, clobber=True)
    test = uvb.read_beamfits(fname) 
    # let's try pickled version
    #    fname=fname.replace(".fits",".pickle")
    #    pickle.dump(uvb,open(fname,'wb'))
    



