#!/bin/bash
uvgen baseunit='-51.0204' stokes='i' systemp='0' harange='-6, 6, 0.1' jyperk='12.7' pnoise='0, 0, 0, 0' telescop='atca' gnoise='0' ant="/Users/Jonathan/Documents/miriad/cat/3.0c.ant" source="/Users/Jonathan/Documents/miriad/cat/point.source" radec='0, -45' corr='0, 1, 0, 104' lat='-30' freq='0.478630092323' out='point_source.vis' ellim='12' cycle='0.1, 0'
uvplt log='uu_vv.log' vis='point_source.vis' device='uuvv.png/PNG' line='wide' options='nobase, log' axis='uu,vv'