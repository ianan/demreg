;; example of how to use the regularized inversion to get DEMs
;; this time using the CHIANTI model DEMs
;; NEED SDO/AIA SolarSoft installed for aia_bp_make_counts_dem.pro
;; i.e. how to use data2dem_reg.pro
;;
;; MODIFICATION HISTORY:
;; 18-11-2011		 Created 		IGH

; Need to make the response functions?
restore,file='aia_resp.dat'

; ; order of filters to use
filt=[0,1,2,3,4,6]
; print,tresp.channels[filt]

; ; data in dn/s/px
; ; Use CHIANTI DEM
; ; options are flare (/fl), active region (/ar) coronal hole (/ch) and quiet sun (/qs)
dem_mod=aia_bp_read_dem(/ar)
dn_in=fltarr(6)
dn_in[0]=aia_bp_make_counts_dem(tresp.a94,/ar)
dn_in[1]=aia_bp_make_counts_dem(tresp.a131,/ar)
dn_in[2]=aia_bp_make_counts_dem(tresp.a171,/ar)
dn_in[3]=aia_bp_make_counts_dem(tresp.a193,/ar)
dn_in[4]=aia_bp_make_counts_dem(tresp.a211,/ar)
dn_in[5]=aia_bp_make_counts_dem(tresp.a335,/ar)


; let's just assume a 1% error throughout
edn_in=0.01*dn_in

; ; Just the temperature response ordered by filters as dn
; ; default units of DN cm^5/ s/px
TRmatrix=tresp.all[*,filt]
; ; what is the logT binning of the temperature responses
logt=tresp.logte

; ;order of regularization, default is 0th
order=0
; ;control the regularization parameter/chisq of result in DEM space: reg_tweak=1
reg_tweak=1
; ;Use guess solution in final regularization? default is no, guess=0.
guess=0.

; run the regularization
reg=data2dem_reg(logT, TRmatrix, dn_in, edn_in,$
  mint=5.2, maxt=6.8, nt=17, $
  order=order,reg_tweak=reg_tweak, guess=guess, $
  channels=tresp.channels[filt])

; now to plot the results....
linecolors
!p.charsize=1.5
; plot the regularized DEM and both vertical and horizontal errors
window,1,xsize=650,ysize=500,title='Regularized DEM'
!p.multi=0
ploterr,reg.logt,reg.dem,reg.elogt,reg.edem,$
  /nohat,errcolor=9,yrange=max(reg.dem)*[1e-2,1.2],$
  xrange=minmax(reg.logt),xstyle=17,$
  xtitle='log!D10!N T',ytitle='DEM(T) [cm!U-5!N K!U-1!N]'
oplot,dem_mod.logte,10d^dem_mod.logdem,color=3,thick=2

window,2,xsize=500,ysize=500,title='DN-space and Residuals'
!p.multi=[0,1,2]
nf=n_elements(reg.channels)
plot,indgen(nf),reg.data,/ylog,psym=6,$
  xrange=[-1,nf],xtickf='(a1)',xticks=nf+1,ystyle=16,thick=3,$
  ytit='DN s!U-1!N',xtit=' ',/nodata,$
  yrange=[0.9*min(reg.data),1.1*max(reg.data)]
oplot,indgen(nf),reg.data_reg,psym=6,color=5,thick=1
oplot,indgen(nf),reg.data,psym=7,color=2,thick=2
for i=0, nf-1 do oplot, [i,i],reg.data[i]+[-reg.edata[i],reg.edata[i]],thick=5,color=2

maxr=1.1*max(abs(reg.residuals))
plot,indgen(nf),reg.residuals,xrange=[-1,nf],xtickn=[' ',reg.channels,' '],$
  xticks=nf+1,ystyle=17,thick=1,yrange=maxr*[-1,1],psym=6,$
  ytit='Residuals',xtit='Filter'
oplot,[-2,nf],[0,0],lines=1
xyouts,-0.5,.75*maxr,'chisq='+string(reg.chisq,format='(f4.1)'),/data

window,3,xsize=500,ysize=500,title='RK Matrix'
!p.multi=0
loadct,8,/silent
gamma_ct,2.
rn=reg.rk
for i=0,n_elements(reg.logt)-1 do rn[i,*]=rn[i,*]/max(rn[i,*])
image_tv,transpose(rn),reg.logt,reg.logt,ytitle='Temperature Bin'
oplot,[0,20],[0,20]
