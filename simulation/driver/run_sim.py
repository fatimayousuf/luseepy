#!/usr/bin/env python
import lusee
import numpy  as np
import healpy as hp
import pickle
import os,sys
import yaml
from yaml.loader import SafeLoader

class SimDriver(dict):
    def __init__ (self,yaml):
        self.update(yaml)
        self.lmax = self['observation']['lmax'] ## common lmax
        self.root = self['paths']['lusee_drive_dir']
        if self.root[0]=='$':
            self.root = os.environ[self.root[1:]]
        self._parse_sky()
        self._parse_beams()
        

    def _parse_sky(self):
        fname = os.path.join(self.root,self['paths']['sky_dir'],self['sky']['file'])
        print ("Loading sky: ",fname)
        self.sky = lusee.sky.FitsSky (fname, lmax = self.lmax)

    def _parse_beams(self):
        broot = os.path.join(self.root,self['paths']['beam_dir'])
        beams = []
        bd = self['beams']
        for b in self['observation']['beams']:
            print ("Loading beam",b,":")
            cbeam = bd[b]
            fname = os.path.join(broot,cbeam['file'])
            print ("  loading file: ",fname)
            B = lusee.LBeam (fname)
            angle = self['observation']['common_beam_angle']+cbeam['angle']
            print ("  rotating: ",angle)
            B.rotate(angle)
            beams.append(B)
        self.beams = beams
        self.Nbeams = len(self.beams)


    def run(self):
        print ("Starting simulation :")
        od = self['observation']
        dt = od['dt']
        if type(dt)==str:
            dt = eval(dt)
        O=lusee.LObservation(od['lunar_day'],deltaT_sec=dt,
                    lun_lat_deg=od['lat'], lun_long_deg = od['long'])
        freq = np.arange(od['freq']['start'],od['freq']['end'],od['freq']['step'])
        print ("  locating frequency indices...")
        freq_list = list(self.sky.freq_list)
        freq_ndx = []
        for f in freq:
            try:
                ndx = freq_list.index(f)
            except ValueError:
                print ("Error:")
                print (f"Frequency {f} does not exist in sky list of frequencies ({freq_list})")
                sys.exit(1)
            freq_ndx.append(ndx)
        print ("  setting up combinations...")
        combs = od['combinations']
        if type(combs)==str:
            if combs=='all':
                combs = []
                for i in range(self.Nbeams):
                    for j in range(i,self.Nbeams):
                        combs.append((i,j))
    
        print ("  setting up Simulation object...")
        S = lusee.Simulator (O,self.beams, self.sky, freq_ndx=freq_ndx, lmax = self.lmax,
                             combinations=combs, Tground = od['Tground'] )
        print ("  Simulating...")
        S.simulate(times=O.times)
        fname = self['simulation']['output']
        print ("Writing to",fname)
        S.write(fname)



if __name__ == "__main__":
    if len(sys.argv)<2:
        print ("Specify yaml config file command line parameter.")
        sys.exit(0)
    yaml_file = sys.argv[1]
    with open(yaml_file) as f:
        config = yaml.load(f,Loader=SafeLoader)
    S=SimDriver(config)
    S.run()
    