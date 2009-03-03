SFR_read_me.txt

Modifications to the Streamflow Routing Package originally documented
by Prudic and others (2004) and Niswonger and Prudic (2005):


Update to the SFR7 Package for MODFLOW-2005 version 1.6.01 (March 2009):

    Following the modification of SFR2 to simulate transient  
streamflow routing based on the kinematic-wave equation, three 
new input variables were added if the option for transient 
routing is used. The option for simulating transient streamflow
routing is documented by Markstrom and others (2008). The following 
three input variables are now required when transient routing
is simulated (IRTFLG>0). NUMTIM, WEIGHT, and FLWTOL follow directly 
after IRTFLG in item 1.

NUMTIM-- is the number of sub time steps used to route 
streamflow. The time step that will be used to route streamflow
will be equal to the MODFLOW time step divided by NUMTIM.

WEIGHT-- is the time weighting factor used to calculate the 
change in channel storage. WEIGHT has a value between 0.5 and 1.
Please refer to equation 83 of Markstrom and others (2008)
for further details.

FLWTOL-- is the streamflow tolerance for convergence of the kinematic
wave equation used for transient streamflow routing. A value of
0.00003 cubic meters per second has been used successfully in test
simulations (and would need to be converted to whatever units are
being used in the particular simulation).


Update to the SFR7 Package for MODFLOW-2005 version 1.6 (January 2009):

    Following the modification of SFR2 to simulate transient  
streamflow routing based on the kinematic-wave equation, two 
new input variables were added if the option for transient 
routing is used. The option for simulating transient streamflow
routing is documented by Markstrom and others (2008). The following 
two input variables are now required when transient routing
is simulated (IRTFLG>0). NUMTIM and WEIGHT follow directly after 
and on the same line as the variable IRTFLG.

NUMTIM-- is the number of sub-time steps used to route 
streamflow. The time step that will be used to route streamflow
will be equal to the MODFLOW time step divided by NUMTIM.

WEIGHT-- is the time weighting factor used to calculate the 
change in channel storage. WEIGHT has values between 0.5 and 1.
Refer to equation 83 of Markstrom and others (2008)
for further details.


Changes to the SFR7 package for MODFLOW-2005 version 1.5:

A few small changes were made to the Streamflow Routing Package
(Niswonger and Prudic, 2005) since its last release. A change was made in
the calculation of stream seepage during steady-state simulations. In the
previous version, single precision variables were used during intermediate
calculations of streambed seepage. As a result, the model could fail to
converge for some steady-state simulations. The code was modified such
that only double precision arithmetic is used during the calculations of
streambed seepage during steady-state simulations.

Another larger modification was made to SFR7 to include the capability of
distributed streamflow routing using the kinematic-wave approximation to
the Saint-Venant Equations (Lighthill and Whitham, 1955). This new
capability is described in Markstrom and others (2008) on pages 68-69.
This capability requires the addition of one new input variable that needs
to be appended to the end of the first record of input variables in the
SFR7. This new input variable is a flag that specifies whether or not the
kinematic-wave equation will be used to route water in channels. The
updated input instructions for SFR7 are documented in Markstrom and others
(2008), pp 202-210.



Update (June 2006):

	Several minor coding changes were made since the initial
release to fix bugs in relation to the different input options for
SFR2 (Niswonger and Prudic, 2005). These fixes to SFR2 resulted
from inconsistencies in the original data input by stream segments
in the SFR1 documentation (Prudic and others, 2004) and the new
data input option by stream reaches in the instructions published by
Niswonger and Prudic (2005). The SFR2 documentation report has
been revised to better explain the different options available in
SFR2. The latest version of the SFR2 documentation report is
version 1.10 and is available in PDF form at the persistent URL
http://pubs.water.usgs.gov/tm6A13/. Corrections to the original
printed document are listed in file tma6a13_SFR2revision_history.pdf
included in the document (doc) directory distributed with MODFLOW.

	An important change was made to SFR2 in the Formulate
and Budget modules that pertain to the computation of outflow from
lakes when streams are connected to lakes in the LAKE(LAK3) Package
(Merritt and Konikow, 2000). These changes were necessary to
remain compatible with the most recent version of the Lake Package.
The initial version of SFR1 and SFR2 made the computation of lake
outflow on the basis of either the lake stage from the previous
time step or the previous MODFLOW iteration or a combination of
both. This formulation of lake outflow can produce an oscillation
in the lake outflow that affects streamflow leakage downstream of
the lake and could prevent MODFLOW from reaching convergence during
a time step. 

	A new subroutine named GWF1SFR2LAKOUTFLW was added to SFR2.
The new subroutine computes the relation of stream stage with
streamflow at the beginning of a stream segment that receives outflow
from a lake. The relation between stream stage and streamflow are
saved in tables that are passed to the revised Lake (LAK3) Package
where the tables are used in computing lake stage and lake outflow to
the stream segment using the Newton iteration method.

	The changes to the Formulate and Budget modules in SFR2
do not affect previous model results unless outflow from a
lake is simulated as inflow to a stream. Model results when
lake outflow is simulated as inflow to a stream could differ
from eariler models that used the previous method of computing
lake outflow. The greatest differences will occur for steady-state
simulations with computed lake outflows to streams or for transient
simulations when a time weighting factor (THETA) of 0.0 was used
for computing lake stage and lake outflow when time steps were long
and lake outflow was sensitive to small changes in lake stage.



References:


Lighthill, M.J., and Whitham, G.B., 1955, On kinematic floods—flood
movements in long rivers: Proceedings, R. Soc. London, v. A220,
p. 281-316.

Markstrom, S.L., Niswonger, R.G., Regan, R.S., Prudic, D.E., and
Barlow, P.M., 2008, GSFLOW—Coupled ground-water and surface-water
flow model based on the integration of the Precipitation-Runoff
Modeling System (PRMS) and the Modular Ground-Water Flow Model
(MODFLOW-2005): U.S. Geological Survey Techniques and Methods 6-D1,
240 p.

Merritt, M.L., and Konikow, L.F., 2000, Documentation of a computer
program to simulate lake-aquifer interaction using the MODFLOW
ground-water model and the MOC3D solute-transport model: U.S.
Geological Survey Water Resources-Investigations Report 00-4167,
146 p.

Niswonger, R.G., and Prudic, D.E.,2005, Documentation of the
Streamflow-Routing (SFR2) Package to include unsaturated flow
beneath streams--A Modification to SFR1: U.S. Geological Survey
Techniques and Methods 6-A13, 48 p.

Prudic, D.E., Konikow, L.F., and Banta, E.R.,2004, A new streamflow-
routing (SFR1) Package to simulate stream-aquifer interaction with
MODFLOW-2000: U.S. Geological Survey Open-File Report 2004-1042, 95 p.


