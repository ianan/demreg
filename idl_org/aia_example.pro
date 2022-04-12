;; example of how to use the regularized inversion to get DEMs
;; i.e. how to use data2dem_reg.pro
;;
;; MODIFICATION HISTORY:
;; 18-11-2011		 Created 		IGH
;; 12-11-2012    Added plotting of DN per T contribution
;; 08-11-2016    Removed ssw dependent routines - should run with standard IDL
;; 20-05-2019    Minor update: must use the timedepend_date option when getting aia response
;; 28-04-2020    Minor tweak to aia_get_response() - added /eve to avoid prompt


; Need to make the response functions?
;if (file_test('aia_resp.dat') eq 0) then $
;  tresp=aia_get_response(/temperature,/dn,/eve,timedepend_date='01-Jul-2010') & save,file='aia_resp.dat',tresp
restore,file='old_aia_respn.dat'

; ; order of filters to use
;; don't want 304
filt=[0,1,2,3,4,6]
; print,tresp.channels[filt]

; ; data in dn/s/px
; ; the order should be the same as filt above
; ; these values are for a Guassian model DEM
dn_in=[274.2, 166.9, 2842.9, 8496.3, 6941.9, 876.1]

; ; error in dn/s/px
; ; taken as the read + photon noise
edn_in=[14.9, 10.9, 41.1, 65.3, 47.1, 14.1]

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
  mint=5.7, maxt=7.3, nt=33, $
  order=order,reg_tweak=reg_tweak, guess=guess, $
  channels=tresp.channels[filt])

; now to plot the results....
loadct,39,/silent
!p.charsize=1.5
; plot the regularized DEM and both vertical and horizontal errors
window,1,xsize=650,ysize=500,title='Regularized DEM'
!p.multi=0
;ploterr,reg.logt,reg.dem,reg.elogt,reg.edem,$
;	/nohat,errcolor=150,yrange=max(reg.dem)*[1e-3,1.2],$
;	xrange=[min(reg.logt),max(reg.logt)],xstyle=17,/ylog,$
;	xtitle='log!D10!N T',ytitle='DEM(T) [cm!U-5!N K!U-1!N]'

plot,reg.logt,reg.dem,yrange=max(reg.dem)*[1e-3,1.2],$
  xrange=[min(reg.logt),max(reg.logt)],xstyle=17,/ylog,$
  xtitle='log!D10!N T',ytitle='DEM(T) [cm!U-5!N K!U-1!N]'

for i=0, n_elements(reg.dem)-1 do oplot,reg.logt[i]+(reg.elogt[i]*[-1,1]),reg.dem[i]*[1,1]
for i=0, n_elements(reg.dem)-1 do oplot,reg.logt[i]*[1,1],reg.dem[i]+(reg.edem[i]*[-1,1])


window,2,xsize=500,ysize=500,title='DN-space and Residuals'
!p.multi=[0,1,2]
loadct,39,/silent
nf=n_elements(reg.channels)
plot,indgen(nf),reg.data,/ylog,psym=6,$
  xrange=[-1,nf],xtickf='(a1)',xticks=nf+1,ystyle=16,thick=3,$
  ytit='DN s!U-1!N',xtit=' ',/nodata,$
  yrange=[0.9*min(reg.data),1.1*max(reg.data)]
oplot,indgen(nf),reg.data_reg,psym=6,color=150,thick=1
;oplot,indgen(nf),reg.data,psym=7,thick=2
for i=0, nf-1 do oplot, [i,i],reg.data[i]+(reg.edata[i]*[-1,1]),thick=1

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

window,4,xsize=600,ysize=400,title='Data Contribution per T'
!p.multi=[0,3,2]
loadct,0,/silent
for i=0, nf-1 do plot,reg.logt,reg.data_cont_t[*,i],$
  title=reg.channels[i]+': '+string(reg.data_reg[i]),xtitle=' log!D10!NT [K]',ytitle='DN s!U-1!N px!U-1!N'

