#!/bin/csh

mkdir PlotFiles

set BINDIR = "/Users/bbradley/Documents/Software/Bin"

set MAIN_TITLE = "GP14.3"

set TYPES = ( slip )
set TYPES = ( slip trise )
set TYPES = ( slip trise rake )
set LABS = ( "Slip (cm)" "Rise Time (s)" "Rake (deg)" )
set GMT_CPTDIR = "./gmt_cpt"
set BCPTS = ( {$GMT_CPTDIR}/y2r2b.cpt {$GMT_CPTDIR}/c2b2m.cpt {$GMT_CPTDIR}/rain6.cpt )

set NEARV = ( 100 0.5 5 )
#the three sets of values below [MAXV, MINV, INCV, CPTA] are the max,min, increments for the color bar, and CPTA is the increment for the legend numbering
set MAXV = ( 800 16.0  600 )
#set MAXV = ( 320 5.0  600 )
#set MAXV = ( 200 3.0  600 )
set MINV = ( 0 0 -600 )
set INCV = ( 80 1.6  10. )
#set INCV = ( 32 0.5  10. )
#set INCV = ( 20 0.3  10. )
set CPTA = ( 400 8.0 20 )
#set CPTA = ( 64 1.0 20 )
#set CPTA = ( 40 0.6 20 )
set RMEAN = ( 0 0 1 )
set CPTZ = ( " " " " " " )
set CPTZ = ( "-Z" "-Z" " " )
set TINIT = ( 1 0 0 )
set RAKES = ( 0 0 1 )
set PSH = ( -0.1 -0.1 -0.3 )
set PREC = ( 0 1 0 )
set SLIP_WGT = ( 0 1 0 )
set COLR_BAR = ( 1 1 0 )
# if the variable below = 1 then the rake angles are defined over [-180,180], otherwise if=0 then angles = [0,360].  This affects the min/avg/max values.  The first two values are for slip and rise time, so shouldnt be adjusted from 0, only the third.  This should generally be "1" when the average rake is near zero.
set FLIP_RAKE_ANGLES = (0 0 1)

set XYZCODE = ( srf2xyz srf2xyz srf2xyz )
set DX     = ( 0.5 0.5 3.3333333 2.0 2.0 )
set DY     = ( 0.5 0.5 2.5000000 2.5 2.5 )
#set DX     = ( 0.2 0.1666667 3.3333333 2.0 2.0 )
#set DY     = ( 0.2 0.1666667 2.5000000 2.5 2.5 )
#set DX     = ( 0.1 0.1666667 3.3333333 2.0 2.0 )
#set DY     = ( 0.1 0.1666667 2.5000000 2.5 2.5 )



# Name of the output file
set SOURCES = ( m7.80-411.0x13.9_s1129571 )

set SLIPDIR = ./
#name (and dir) of the gmt ps file
set SLIPS = ( Srf/m7.80-411.0x13.9_s1129571.srf )

#fault geometry - length and downdip width
set FLEN   = ( 411.0 )
set FWID   = ( 13.9 )

set DX_RAKE = 2.0
set DY_RAKE = 2.0
#set DX_RAKE = 0.67
#set DY_RAKE = 0.67
set USE_AVG_RAKE = 1
set USE_AVG_RAKE = 0

#set KMINCH = ( 4.0 )  #for Mw~5
#set KMINCH = ( 8.0 )  #Mw~6
#set KMINCH = ( 12.0 ) #Mw~7
set KMINCH = ( 40.0 ) #Mw~8


set NCOL   = (   6 )
#XTIC YTIC are the axis increment ticks to plot
set XTIC   = (   20 )
set YTIC   = (   6 )
#shift of plot
#set YZ     = ( 1.7 )
set YZ     = ( 5.5 )
#This is the contour interval for the rupt time contours
set CINTV  = ( 2.0 )
#set CINTV  = ( 1.0 )

set XPMID = 4.25
set CALC_XY = 0

set GRDFILE = r.grd

#shift
set XZ = 0.5
set XSHFT = 3.2

set XLAB = "along strike (km)"
set YLAB = "W (km)"


#gmt gmtset FRAME_WIDTH 0.05 LABEL_FONT_SIZE 11 ANOT_FONT_SIZE 11 PLOT_DEGREE_FORMAT D PAGE_ORIENTATION PORTRAIT TICK_LENGTH 0.03 D_FORMAT %lg
gmt gmtset MAP_FRAME_WIDTH 0.05 FONT_LABEL 11 FORMAT_GEO_MAP D MAP_TICK_LENGTH_PRIMARY 0.03 FORMAT_FLOAT_OUT %lg PROJ_LENGTH_UNIT inch PS_PAGE_ORIENTATION LANDSCAPE

set m = 0
foreach source ( $SOURCES )
@ m ++

set PSFILE = `echo $source | gawk '{printf "PlotFiles/%s.ps\n",$1;}'`
\rm $PSFILE

set XSH = 0.0
set YSH = 0.0

# FIND MAX. SLIP FOR ALL SLIPMODELS

set TCPT = ( )
set AINC = ( )
set TINC = ( )
set t = 0
foreach typ ( $TYPES )
@ t ++

set TCPT = ( $TCPT temp${t}.cpt )

set SMIN = 1.0e+15
set SMAX = -1.0e+15
set s = 0
foreach slip ( $SLIPS )
@ s ++

#set SMAX = `{$BINDIR}/$XYZCODE[$s] calc_xy=$CALC_XY type="${typ}" nseg=-1 < $SLIPDIR/${slip} | gawk -v m=$SMAX -v p=$PREC[$t] '{if($3>m)m=$3;}END{fmt=sprintf("%%.%df\n",p);printf fmt,m;}'`
#set SMIN = `{$BINDIR}/$XYZCODE[$s] calc_xy=$CALC_XY type="${typ}" nseg=-1 < $SLIPDIR/${slip} | gawk -v m=$SMIN -v p=$PREC[$t] '{if($3<m)m=$3;}END{fmt=sprintf("%%.%df\n",p);printf fmt,m;}'`
#the two lines below include the possibility of negative rakes (see FLIP_RAKE_ANGLES parameter)
set SMAX = `{$BINDIR}/$XYZCODE[$s] calc_xy=$CALC_XY type="${typ}" nseg=-1 < $SLIPDIR/${slip} | gawk -v m=$SMAX -v p=$PREC[$t] -v fliprake=$FLIP_RAKE_ANGLES[$t] '{val=$3;if(fliprake==1){if(val>180.0)val=val-360.0;}if(val>m)m=val;}END{fmt=sprintf("%%.%df\n",p);printf fmt,m;}'`
set SMIN = `{$BINDIR}/$XYZCODE[$s] calc_xy=$CALC_XY type="${typ}" nseg=-1 < $SLIPDIR/${slip} | gawk -v m=$SMIN -v p=$PREC[$t] -v fliprake=$FLIP_RAKE_ANGLES[$t] '{val=$3;if(fliprake==1){if(val>180.0)val=val-360.0;}if(val<m)m=val;}END{fmt=sprintf("%%.%df\n",p);printf fmt,m;}'`

end

echo $SMIN $SMAX

set MAXRND = `echo $SMAX | gawk -v p=$PREC[$t] -v nv=$NEARV[$t] '{m=$1*1.1;}END{d=int($1/nv+0.0);while(nv*d<0.9*m)d++;fmt=sprintf("%%.%df %%.2f\n",p);printf fmt,nv*d,0.2*nv*d;}'`
#echo $MAXRND

if($MAXV[$t] != "-1") then

set MAXRND = ( $MAXV[$t] $INCV[$t] )
set SMIN = $MINV[$t]

set AA = $CPTA[$t]

else

set AA = ${MAXRND[2]}

endif

set AINC = ( $AINC $AA )
set TINC = ( $TINC ${MAXRND[2]} )

gmt makecpt -C$BCPTS[$t] -T${SMIN}/$MAXRND[1]/$MAXRND[2] $CPTZ[$t] > $TCPT[$t]

end

#removed from pstext "-G0/0/0"
gmt pstext -JX8.5/11.0 -R0/8.5/0/11 -N -K -X0.0 -Y0.0 << END > $PSFILE
#1.2 1.0 20 0 1 1 ${source}
#1.2 0.8 10 0 0 1 ${SLIPS[1]}
END

echo $cwd
#set LAB = "-U/0.2/0.2/$cwd"
#
#pstext $LAB -JX8.5/11.0 -R0/8.5/0/11 -N -G0/0/0 -O -K << END >> $PSFILE
#4.25 9.80 12 0 1 6 Northridge Slipmodels
#END

#removed from pstext "-G0/0/0"
gmt pstext -JX8.5/11.0 -R0/8.5/0/11 -N -O -K -X$XZ -Y$YZ[$m] << END >> $PSFILE
END

set s = 0
foreach slip ( $SLIPS )
@ s ++

set XINCH = `echo $FLEN[$s] ${KMINCH[$m]} | gawk '{printf "%f\n",$1/$2;}'`
set YINCH = `echo $FWID[$s] ${KMINCH[$m]} | gawk '{printf "%f\n",$1/$2;}'`
set XMID = `echo $XINCH | gawk '{printf "%f\n",$1+0.3;}'`
set YSLEN = `echo $YINCH | gawk '{printf "%f\n",0.8*$1;}'`
set YSMID = `echo $YSLEN | gawk '{printf "%f\n",0.5*$1;}'`

set SCALE = "$XINCH/-$YINCH"

set REGION = "0/$FLEN[$s]/0/$FWID[$s]"
set ATTRIB = "-JX$SCALE -R$REGION"

set t = 0
foreach typ ( $TYPES )
@ t ++

if( $t == 1 ) then

#removed from pstext "-G0/0/0"
echo $FLEN[$s] | gawk -v t="${MAIN_TITLE}" '{printf "%f 0 20 0 1 2 %s\n",0.5*$1,t;}' | \
gmt pstext -N $ATTRIB -D0.0/0.2  -K -O >> $PSFILE

endif

set W = "w"
if( $s <= ${NCOL[$m]} ) then
   set W = "W"
endif

set S = "s"
if( ($s % ${NCOL[$m]}) == 0 || (($s == $#SLIPS) && ($t == $#TYPES))) then
   set S = "S"
endif

set SLIPFILE = $SLIPDIR/${slip}

#set AVG_MAX = `{$BINDIR}/$XYZCODE[$s] calc_xy=$CALC_XY type="${typ}" nseg=-1 < $SLIPFILE | gawk -v p=$PREC[$t] 'BEGIN{mx=-1.0e+15;mn=1.0e+15;}{s=s+$3;if($3>mx)mx=$3;if($3<mn)mn=$3;}END{fmt=sprintf("%%.%df %%.%df %%.%df\n",p,p,p);printf fmt,s/NR,mx,mn;}'`
#set AVG_MAX = `{$BINDIR}/$XYZCODE[$s] calc_xy=$CALC_XY type="${typ}" nseg=-1 dump_slip=1 < $SLIPFILE | gawk -v p=$PREC[$t] -v sw=$SLIP_WGT[$t] 'BEGIN{mx=-1.0e+15;mn=1.0e+15;}{w=1;if(sw==1)w=$4;v=v+$3*w;tw=tw+w;if($3>mx)mx=$3;if($3<mn)mn=$3;}END{fmt=sprintf("%%.%df %%.%df %%.%df\n",p,p,p);printf fmt,v/tw,mx,mn;}'`
set AVG_MAX = `{$BINDIR}/$XYZCODE[$s] calc_xy=$CALC_XY type="${typ}" nseg=-1 dump_slip=1 < $SLIPFILE | gawk -v p=$PREC[$t] -v sw=$SLIP_WGT[$t] -v fliprake=$FLIP_RAKE_ANGLES[$t] 'BEGIN{mx=-1.0e+15;mn=1.0e+15;}{w=1;if(sw==1)w=$4;val=$3;if(fliprake==1){if(val>180.0)val=val-360.0;}v=v+val*w;tw=tw+w;if(val>mx)mx=val;if(val<mn)mn=val;}END{fmt=sprintf("%%.%df %%.%df %%.%df\n",p,p,p);printf fmt,v/tw,mx,mn;}'`
echo $s avg= $AVG_MAX[1] max= $AVG_MAX[2] min= $AVG_MAX[3]

#changed -F in 'xyz2grd to -r
{$BINDIR}/$XYZCODE[$s] calc_xy=$CALC_XY type="${typ}" nseg=-1 < $SLIPFILE | \
gawk -v rmean=$RMEAN[$t] -v avg=$AVG_MAX[1] -v fliprake=$FLIP_RAKE_ANGLES[$t] 'BEGIN{vv=0.0;if(rmean==1)vv=avg;}{ \
val=$3;if(fliprake==1){if(val>180.0)val=val-360.0;} \
printf "%13.5e %13.5e %13.5e\n",$1,$2,val-vv;}' | \
gmt xyz2grd -G$GRDFILE -I$DX[$s]/$DY[$s] -R$REGION -r
#plot the grid image
gmt grdimage $GRDFILE $ATTRIB -C$TCPT[$t] -B${XTIC[$m]}:"$XLAB":/${YTIC}:"$YLAB":${W}${S}en -K -O -X$XSH -Y$YSH >> $PSFILE

if($TINIT[$t] == 1) then

#changed -F in 'xyz2grd to -r
{$BINDIR}/$XYZCODE[$s] calc_xy=$CALC_XY type="tinit" nseg=-1 < $SLIPFILE | \
gmt xyz2grd -G$GRDFILE -I$DX[$s]/$DY[$s] -R$REGION -r
gmt grdcontour $GRDFILE $ATTRIB -C$CINTV[$m] -W1.0  -K -O >> $PSFILE
#removed '-W4/0/0/0' from the above line

endif

if($RAKES[$t] == 1) then

set SMAX = `{$BINDIR}/$XYZCODE[$s] calc_xy=$CALC_XY type=slip nseg=-1 < $SLIPFILE | gawk -v p=$PREC[$t] '{if($3>m)m=$3;}END{fmt=sprintf("%%.%df\n",p);printf fmt,m;}'`

#because YDIR is reversed then 'rk[i]/nv[i]' changed to '-rk[i]/nv[i]' (near end of pipe below)
{$BINDIR}/$XYZCODE[$s] calc_xy=$CALC_XY type=rake nseg=-1 dump_slip=1 < $SLIPFILE | \
gawk -v dx=$DX_RAKE -v dy=$DY_RAKE -v len=$FLEN[$s] -v wid=$FWID[$s] -v avgr=$USE_AVG_RAKE -v mx=$SMAX 'BEGIN{ \
nx=int(len/dx+0.5);ny=int(wid/dy+0.5);for(i=1;i<=nx*ny;i++){mr[i]=1.0e+15;x0[i]=0.0;y0[i]=0.0;nv[i]=0;rk[i]=0.0;sp[i]=0.0;}}{ \
ix=int($1/dx);iy=int($2/dy);ip=1+ix+iy*nx; \
if(avgr==0){ \
xx=(ix+0.5)*dx - $1; yy=(iy+0.5)*dy - $2;if((xx*xx+yy*yy)<mr[ip]) { \
x0[ip]=(ix+0.5)*dx; y0[ip]=(iy+0.5)*dy;nv[ip]=1;rk[ip]=$3;mr[ip]=(xx*xx+yy*yy);sp[ip]=$4;}} \
else{ \
x0[ip]=(ix+0.5)*dx; y0[ip]=(iy+0.5)*dy;nv[ip]++;rk[ip]=rk[ip]+$3;sp[ip]=sp[ip]+$4;} \
} \
END{ \
for(i=1;i<=nx*ny;i++){ \
if(nv[i]>0)printf "%13.5e %13.5e %13.5e %f\n",x0[i],y0[i],-rk[i]/nv[i],0.4*sp[i]/(nv[i]*mx);}}' | \
gmt psxy $ATTRIB -Sv0.005i/0.04i/0.02i -W2 -G0/0/0 -K -O >> $PSFILE
# changed the line above from that below:
# gmt psxy $ATTRIB -Sv0.005i/0.04i/0.02i -W1/0/0/0 -G0/0/0 -K -O >> $PSFILE
endif

#removed from pstext "-G0/0/0"
gmt pstext $ATTRIB -N -O -K -D0.025/0.05 << END >> $PSFILE
#0.0 0.0 11 0 0 1 $slip \b\b\b $AVG_MAX[1] / $AVG_MAX[2]
0.0 0.0 12 0 1 1 $LABS[$t]
$FLEN[$s] 0.0 11 0 0 3 $AVG_MAX[3] / $AVG_MAX[1] / $AVG_MAX[2]
END

if( $COLR_BAR[$t] == 1 ) then

gmt psscale -Ef -C$TCPT[$t] -D${XMID}/${YSMID}/${YSLEN}/0.1 -K -O -Ba${AINC[$t]}f${TINC[$t]}g${TINC[$t]}::/:: >> $PSFILE

endif

set XSH = 0.0
if( ($s % ${NCOL[$m]}) == 0) then
   set XSH = `echo $XINCH | gawk '{printf "%f\n",$1+0.2;}'`
endif

set YSH = `echo $YINCH | gawk '{printf "%f\n",-($1+0.4);}'`
if( ($s % ${NCOL[$m]}) == 0) then
   set YSH = `echo $YSH ${NCOL[$m]} | gawk '{printf "%f\n",-$1*($2-1);}'`
endif

end

end

gmt psxy $ATTRIB -W5 -G0/0/0 -O << END >> $PSFILE
END

end

\rm $GRDFILE $TCPT
