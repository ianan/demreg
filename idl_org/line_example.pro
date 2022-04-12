;; example of how to use the regularized inversion to get DEMs
;; i.e. how to use data2dem_reg.pro
;;
;; MODIFICATION HISTORY:
;; 18-11-2011		 Created 		IGH
;; 08-11-2016    Removed ssw dependent routines - should run with standard IDL

; ; Load some contribution functions
;; here are just 9 lines taken from Landi & Young 2009
;; if you want your own try using gofnt.pro
;restgen,file='con_func_9',cf
restore,file='con_func_9.dat'

; ; line intensities in erg/cm^2/s/sr
;; these intesities for the lines in cf are from the CHIANTI /qs DEM model
line_in=[489.7,  267.5, 4.9,  77.4, 67.2,  26.5,  47.1,  235.9, 6.9]

; ; error in erg/cm^2/s/sr
; ; just assume % of line intensities
eline_in=line_in*0.01

; ; get all the contribution functions from the save structure
; ;  units of are erg cm^3/s/sr
CFmatrix=cf.con_func
; ; what is the logT binning of the contribution functions
logt=cf[0].logt

; ;order of regularization, default is 0th
order=0
; ;control the regularization parameter/chisq of result in DEM space: reg_tweak=1
reg_tweak=1
; ;Use guess solution in final regularization? default is no, guess=0.
guess=0
;; Use the min of the EM loci curves as the initial guess solution
;; used to weight/create the constraint matrix and possibly in the regularization itself (if guess=1)
gloci=0
;; Also choose reg parameter to that gives a positive DEM?
pos=1

; run the regularization
reg=data2dem_reg(logT, CFmatrix, line_in, eline_in,$
  mint=5.0, maxt=6.7, nt=50, $
  order=order,reg_tweak=reg_tweak, guess=guess, $
  channels=cf.line,gloci=gloci,pos=pos)

; now to plot the results....
loadct,39
!p.charsize=1.5
; plot the regularized DEM and both vertical and horizontal errors
window,1,xsize=650,ysize=500,title='Regularized DEM'
!p.multi=0
;ploterr,reg.logt,reg.dem,reg.elogt,reg.edem,$
; /nohat,errcolor=150,yrange=max(reg.dem)*[1e-3,1.2],$
; xrange=[min(reg.logt),max(reg.logt)],xstyle=17,/ylog,$
; xtitle='log!D10!N T',ytitle='DEM(T) [cm!U-5!N K!U-1!N]'

plot,reg.logt,reg.dem,yrange=max(reg.dem)*[1e-3,1.2],$
  xrange=[min(reg.logt),max(reg.logt)],xstyle=17,/ylog,$
  xtitle='log!D10!N T',ytitle='DEM(T) [cm!U-5!N K!U-1!N]'

for i=0, n_elements(reg.dem)-1 do oplot,reg.logt[i]+(reg.elogt[i]*[-1,1]),reg.dem[i]*[1,1]
for i=0, n_elements(reg.dem)-1 do oplot,reg.logt[i]*[1,1],reg.dem[i]+(reg.edem[i]*[-1,1])
;dem_mod=aia_bp_read_dem(/qs)
;oplot,dem_mod.logte,10d^dem_mod.logdem,color=3,thick=2

window,2,xsize=500,ysize=500,title='Line Intensities and Residuals'
!p.multi=[0,1,2]
loadct,39,/silent
nf=n_elements(reg.channels)
plot,indgen(nf),reg.data,/ylog,psym=6,$
  xrange=[-1,nf],xtickf='(a1)',xticks=nf+1,ystyle=16,thick=3,$
  ytit='Line Intensities',xtit=' ',/nodata,$
  yrange=[0.9*min(reg.data),1.1*max(reg.data)]
oplot,indgen(nf),reg.data_reg,psym=6,color=150,thick=1
;oplot,indgen(nf),reg.data,psym=7,thick=2
for i=0, nf-1 do oplot, [i,i],reg.data[i]+(reg.edata[i]*[-1,1]),thick=1

maxr=1.1*max(abs(reg.residuals))
plot,indgen(nf),reg.residuals,xrange=[-1,nf],xtickn=[' ',reg.channels,' '],$
  xticks=nf+1,ystyle=17,thick=1,yrange=maxr*[-1,1],psym=6,$
  ytit='Residuals',xtit='Line'
oplot,[-2,nf],[0,0],lines=1
xyouts,-0.5,.75*maxr,'chisq='+string(reg.chisq,format='(f4.1)'),/data

window,3,xsize=800,ysize=500,title='RK Matrix'
!p.multi=[0,2,1]
loadct,8,/silent
gamma_ct,2.
rn=reg.rk
for i=0,n_elements(reg.logt)-1 do rn[i,*]=rn[i,*]/max(rn[i,*])
image_tv,transpose(rn),reg.logt,reg.logt,ytitle='Temperature Bin'
loadct,39,/silent
oplot,[0,20],[0,20]
oplot,[0,20],reg.logt[25]*[1.,1.],lines=2,color=100,thick=3
oplot,[0,20],reg.logt[11]*[1.,1.],lines=1,color=150,thick=3
oplot,[0,20],reg.logt[38]*[1.,1.],lines=3,color=200,thick=3
plot,reg.logt,rn[25,*],xstyle=17,/nodata,yrange=[0,1.1]
oplot,reg.logt,rn[25,*],color=100,lines=2,thick=3
oplot,reg.logt,rn[11,*],color=150,lines=1,thick=3
oplot,reg.logt,rn[38,*],color=200,lines=3,thick=3