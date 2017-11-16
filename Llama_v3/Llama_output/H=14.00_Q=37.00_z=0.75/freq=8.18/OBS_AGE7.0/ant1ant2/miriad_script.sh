#!/bin/bash
uvgen baseunit='-51.0204' stokes='i' systemp='0' harange='-6, 6, 0.1' jyperk='12.7' pnoise='0, 0, 0, 0' telescop='atca' gnoise='0' ant="/Users/Jonathan/Documents/miriad/cat/3.0c.ant" source="/Users/Jonathan/Documents/miriad/cat/point.source" radec='0, -45' corr='0, 1, 0, 104' lat='-30' freq='0.151356124844' out='point_source.vis' ellim='12' cycle='0.1, 0'
fits op='xyin' out='raise.mir' in='../../age=7.0_z=0.75_freq=8.18_beam_d=1.0.fits'
uvmodel vis='point_source.vis' model='raise.mir' options='replace' out='raise.vis'
uvplt log='uvdistance_phase.log' vis='raise.vis' device='uvdistancephase.png/PNG' line='wide' options='log' select='ant(1)(2)' axis='uvdistance,phase'
uvplt log='uvdistance_amplitude.log' vis='raise.vis' device='uvdistanceamplitude.png/PNG' line='wide' options='log' select='ant(1)(2)' axis='uvdistance,amplitude'
uvplt log='uvdistance_dtime.log' vis='raise.vis' device='uvdistancedtime.png/PNG' line='wide' options='log' select='ant(1)(2)' axis='uvdistance,dtime'