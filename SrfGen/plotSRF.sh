#!/usr/bin/env bash

mkdir -p PlotFiles

srf2xyz_bin="/home/vap30/bin/srf2xyz"

title="GP14.3"

TYPES=( slip trise rake )
LABS=( "Slip (cm)" "Rise Time (s)" "Rake (deg)" )
cpt_dir="./srf_cpt"
BCPTS=( "$cpt_dir/y2r2b.cpt" "$cpt_dir/c2b2m.cpt" "$cpt_dir/rain6.cpt" )

NEARV=( 100 0.5 5 )
#the three sets of values below [MAXV, MINV, INCV, CPTA] are the max,min, increments for the color bar, and CPTA is the increment for the legend numbering
MAXV=( 800 16.0  600 )
MINV=( 0 0 -600 )
INCV=( 80 1.6  10. )
CPTA=( 400 8.0 20 )
RMEAN=( 0 0 1 )
#CPTZ = ( " " " " " " )
CPTZ=( "-Z" "-Z" " " )
TINIT=( 1 0 0 )
RAKES=( 0 0 1 )
PSH=( -0.1 -0.1 -0.3 )
PREC=( 0 1 0 )
SLIP_WGT=( 0 1 0 )
COLR_BAR=( 1 1 0 )
# if the variable below = 1 then the rake angles are defined over [-180,180], otherwise if=0 then angles = [0,360].  This affects the min/avg/max values.  The first two values are for slip and rise time, so shouldnt be adjusted from 0, only the third.  This should generally be "1" when the average rake is near zero.
FLIP_RAKE_ANGLES=(0 0 1)

DX=( 0.5 0.5 3.3333333 2.0 2.0 )
DY=( 0.5 0.5 2.5000000 2.5 2.5 )

# Name of the output file
SOURCES=( bev01_s103245 )

SLIPDIR=./
#name (and dir) of the gmt ps file
SLIPS=( Srf/bev01_s103245.srf )

#fault geometry - length and downdip width
FLEN=( 411.0 )
FWID=( 13.9 )

DX_RAKE=2.0
DY_RAKE=2.0
USE_AVG_RAKE=0

#KMINCH = ( 4.0 )  #for Mw~5
#KMINCH = ( 8.0 )  #Mw~6
#KMINCH = ( 12.0 ) #Mw~7
KMINCH=( 40.0 ) #Mw~8


NCOL=(   6 )
#XTIC YTIC are the axis increment ticks to plot
XTIC=(   20 )
YTIC=(   6 )
#shift of plot
#YZ     = ( 1.7 )
YZ=( 5.5 )
#This is the contour interval for the rupt time contours
CINTV=( 2.0 )
#CINTV  = ( 1.0 )

XPMID=4.25
CALC_XY=0

GRDFILE=r.grd

#shift
XZ=0.5
XSHFT=3.2

XLAB="along strike (km)"
YLAB="W (km)"


gmtset MAP_FRAME_WIDTH 0.05 FONT_LABEL 11 FORMAT_GEO_MAP D MAP_TICK_LENGTH_PRIMARY 0.03 FORMAT_FLOAT_OUT %lg PROJ_LENGTH_UNIT inch PS_PAGE_ORIENTATION LANDSCAPE

m=0
for source in $SOURCES; do
    PSFILE=`echo $source | gawk '{printf "PlotFiles/%s.ps\n",$1;}'`

    \rm $PSFILE

    XSH=0.0
    YSH=0.0

    # FIND MAX. SLIP FOR ALL SLIPMODELS
    TCPT=( )
    AINC=( )
    TINC=( )
    t=0
    for typ in $TYPES; do

        TCPT=( $TCPT temp${t}.cpt )

        SMIN=1.0e+15
        SMAX=-1.0e+15
        s=0
        for slip in $SLIPS; do
            #the two lines below include the possibility of negative rakes (see FLIP_RAKE_ANGLES parameter)
            SMAX=`$srf2xyz_bin calc_xy=$CALC_XY type="${typ}" nseg=-1 < $SLIPDIR/${slip} | gawk -v m=$SMAX -v p=${PREC[$t]} -v fliprake=$FLIP_RAKE_ANGLES[$t] '{val=$3;if(fliprake==1){if(val>180.0)val=val-360.0;}if(val>m)m=val;}END{fmt=sprintf("%%.%df\n",p);printf fmt,m;}'`
            SMIN=`$srf2xyz_bin calc_xy=$CALC_XY type="${typ}" nseg=-1 < $SLIPDIR/${slip} | gawk -v m=$SMIN -v p=${PREC[$t]} -v fliprake=${FLIP_RAKE_ANGLES[$t]} '{val=$3;if(fliprake==1){if(val>180.0)val=val-360.0;}if(val<m)m=val;}END{fmt=sprintf("%%.%df\n",p);printf fmt,m;}'`
            s=$(($s + 1))
        done

        echo $SMIN $SMAX

        MAXRND=`echo $SMAX | gawk -v p=${PREC[$t]} -v nv=${NEARV[$t]} '{m=$1*1.1;}END{d=int($1/nv+0.0);while(nv*d<0.9*m)d++;fmt=sprintf("%%.%df %%.2f\n",p);printf fmt,nv*d,0.2*nv*d;}'`
        #echo $MAXRND

        if [ ${MAXV[$t]} != "-1" ]; then
            MAXRND=( ${MAXV[$t]} {$INCV[$t]} )
            SMIN=${MINV[$t]}
            AA=${CPTA[$t]}
        else
            AA=${MAXRND[1]}
        fi

        AINC=( $AINC $AA )
        TINC=( $TINC ${MAXRND[1]} )

        makecpt -C${BCPTS[$t]} -T${SMIN}/${MAXRND[0]}/{$MAXRND[1]} ${CPTZ[$t]} > $TCPT[$t]
        t=$(($t + 1))
    done

    # removed from pstext "-G0/0/0"
    pstext -JX8.5/11.0 -R0/8.5/0/11 -N -K -X0.0 -Y0.0 << END > $PSFILE
#1.2 1.0 20 0 1 1 ${source}
#1.2 0.8 10 0 0 1 ${SLIPS[0]}
END

    echo $cwd
    #LAB = "-U/0.2/0.2/$cwd"
    #
    #pstext $LAB -JX8.5/11.0 -R0/8.5/0/11 -N -G0/0/0 -O -K << END >> $PSFILE
#4.25 9.80 12 0 1 6 Northridge Slipmodels
#END

    #removed from pstext "-G0/0/0"
    gmt pstext -JX8.5/11.0 -R0/8.5/0/11 -N -O -K -X$XZ -Y${YZ[$m]} << END >> $PSFILE
END

    s=0
    for slip in $SLIPS; do
        XINCH=`echo ${FLEN[$s]} ${KMINCH[$m]} | gawk '{printf "%f\n",$1/$2;}'`
        YINCH=`echo ${FWID[$s]} ${KMINCH[$m]} | gawk '{printf "%f\n",$1/$2;}'`
        XMID=`echo $XINCH | gawk '{printf "%f\n",$1+0.3;}'`
        YSLEN=`echo $YINCH | gawk '{printf "%f\n",0.8*$1;}'`
        YSMID=`echo $YSLEN | gawk '{printf "%f\n",0.5*$1;}'`

        SCALE="$XINCH/-$YINCH"

        REGION="0/${FLEN[$s]}/0/${FWID[$s]}"
        ATTRIB="-JX$SCALE -R$REGION"

        t=0
        for typ in $TYPES; do

            if [ $t -eq 1 ]; then
                #removed from pstext "-G0/0/0"
                echo $FLEN[$s] | gawk -v t="$title" '{printf "%f 0 20 0 1 2 %s\n",0.5*$1,t;}' | \
                pstext -N $ATTRIB -D0.0/0.2  -K -O >> $PSFILE
            fi

            W="w"
            if [ "$s" -le ${NCOL[$m]} ]; then
               W="W"
            fi

            S="s"
            if [ $(($s % ${NCOL[$m]})) -eq 0 ] || ( [ $s == ${#SLIPS[@]} ] && [ $t == ${#TYPES[@]} ] ); then
                S="S"
            fi

            SLIPFILE=$SLIPDIR/${slip}
exit
            AVG_MAX=`$srf2xyz_bin calc_xy=$CALC_XY type="${typ}" nseg=-1 dump_slip=1 < $SLIPFILE | gawk -v p=${PREC[$t]} -v sw=${SLIP_WGT[$t]} -v fliprake=${FLIP_RAKE_ANGLES[$t]} 'BEGIN{mx=-1.0e+15;mn=1.0e+15;}{w=1;if(sw==1)w=$4;val=$3;if(fliprake==1){if(val>180.0)val=val-360.0;}v=v+val*w;tw=tw+w;if(val>mx)mx=val;if(val<mn)mn=val;}END{fmt=sprintf("%%.%df %%.%df %%.%df\n",p,p,p);printf fmt,v/tw,mx,mn;}'`
            echo $s avg= ${AVG_MAX[0]} max= ${AVG_MAX[1]} min= ${AVG_MAX[2]}

            #changed -F in 'xyz2grd to -r
            $srf2xyz_bin calc_xy=$CALC_XY type="${typ}" nseg=-1 < $SLIPFILE | \
            gawk -v rmean=${RMEAN[$t]} -v avg=${AVG_MAX[0]} -v fliprake=$FLIP_RAKE_ANGLES[$t] 'BEGIN{vv=0.0;if(rmean==1)vv=avg;}{ \
                    val=$3;if(fliprake==1){if(val>180.0)val=val-360.0;} \
                    printf "%13.5e %13.5e %13.5e\n",$1,$2,val-vv;}' | \
            xyz2grd -G$GRDFILE -I${DX[$s]}/${DY[$s]} -R$REGION -r
            #plot the grid image
            grdimage $GRDFILE $ATTRIB -C${TCPT[$t]} -B${XTIC[$m]}:"$XLAB":/${YTIC}:"$YLAB":${W}${S}en -K -O -X$XSH -Y$YSH >> $PSFILE

            if [ ${TINIT[$t]} -eq 1 ]; then
                #changed -F in 'xyz2grd to -r
                $srf2xyz_bin calc_xy=$CALC_XY type="tinit" nseg=-1 < $SLIPFILE | \
                xyz2grd -G$GRDFILE -I${DX[$s]}/${DY[$s]} -R$REGION -r
                grdcontour $GRDFILE $ATTRIB -C${CINTV[$m]} -W1.0  -K -O >> $PSFILE
                #removed '-W4/0/0/0' from the above line
            fi

            if [ $RAKES[$t] -eq 1 ]; then
                SMAX = `$srf2xyz_bin calc_xy=$CALC_XY type=slip nseg=-1 < $SLIPFILE | gawk -v p=${PREC[$t]} '{if($3>m)m=$3;}END{fmt=sprintf("%%.%df\n",p);printf fmt,m;}'`

                #because YDIR is reversed then 'rk[i]/nv[i]' changed to '-rk[i]/nv[i]' (near end of pipe below)
                $srf2xyz_bin calc_xy=$CALC_XY type=rake nseg=-1 dump_slip=1 < $SLIPFILE | \
                gawk -v dx=$DX_RAKE -v dy=$DY_RAKE -v len=${FLEN[$s]} -v wid=${FWID[$s]} -v avgr=$USE_AVG_RAKE -v mx=$SMAX 'BEGIN{ \
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
                psxy $ATTRIB -Sv0.005i/0.04i/0.02i -W2 -G0/0/0 -K -O >> $PSFILE
                # changed the line above from that below:
                # gmt psxy $ATTRIB -Sv0.005i/0.04i/0.02i -W1/0/0/0 -G0/0/0 -K -O >> $PSFILE
            fi

            #removed from pstext "-G0/0/0"
            pstext $ATTRIB -N -O -K -D0.025/0.05 << END >> $PSFILE
#0.0 0.0 11 0 0 1 $slip \b\b\b ${AVG_MAX[0]} / ${AVG_MAX[1]}
0.0 0.0 12 0 1 1 ${LABS[$t]}
${FLEN[$s]} 0.0 11 0 0 3 ${AVG_MAX[2]} / ${AVG_MAX[0]} / ${AVG_MAX[1]}
END

            if [ $COLR_BAR[$t] -eq 1 ]; then
                gmt psscale -Ef -C$TCPT[$t] -D${XMID}/${YSMID}/${YSLEN}/0.1 -K -O -Ba${AINC[$t]}f${TINC[$t]}g${TINC[$t]}::/:: >> $PSFILE
            fi

            XSH=0.0
            if [ $(($s % ${NCOL[$m]})) -eq 0 ]; then
                XSH=`echo $XINCH | gawk '{printf "%f\n",$1+0.2;}'`
            fi

            YSH=`echo $YINCH | gawk '{printf "%f\n",-($1+0.4);}'`
            if [ $(($s % ${NCOL[$m]})) -eq 0]; then
                YSH=`echo $YSH ${NCOL[$m]} | gawk '{printf "%f\n",-$1*($2-1);}'`
            fi

            t=$((t+1))
        done
        s=$(($s + 1))
    done

    psxy $ATTRIB -W5 -G0/0/0 -O << END >> $PSFILE
END

    m=$(($m + 1))
done

\rm $GRDFILE $TCPT
