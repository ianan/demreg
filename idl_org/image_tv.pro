pro image_tv,a,xx,yy ,WINDOW_SCALE = window_scale, ASPECT = aspect, $
	INTERP = interp,ytitle=ytitle

; modified version of image_cont.pro just to plot a scaled image IGH



on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if sz[0] lt 2 then message, 'Parameter not 2D'

	;set window used by contour
contour,a,xx,yy,/nodata, xstyle=17, ystyle = 17,ytitle=ytitle,xtitle=xtitle
px = !x.window * !d.x_vsize	;Get size of window in device units
py = !y.window * !d.y_vsize
swx = px[1]-px[0]		;Size in x in device units
swy = py[1]-py[0]		;Size in Y
six = float(sz[1])		;Image sizes
siy = float(sz[2])
aspi = six / siy		;Image aspect ratio
aspw = swx / swy		;Window aspect ratio
f = aspi / aspw			;Ratio of aspect ratios

dosub = 0
if (!d.flags and 1) ne 0 then begin	;Scalable pixels?
  if keyword_set(aspect) then begin	;Retain aspect ratio?
				;Adjust window size
	if f ge 1.0 then swy = swy / f else swx = swx * f
	endif

  tvscl,a,px[0],py[0],xsize = swx, ysize = swy, /device

endif else begin	;Not scalable pixels	
   if keyword_set(window_scale) then begin ;Scale window to image?
	tvscl,a,px[0],py[0]	;Output image
	swx = six		;Set window size from image
	swy = siy
        ; If image is larger than window, only contour portion
        ; that fits within the window.
        doSubx = (px[0]+swx) gt !d.x_size
        doSuby = (py[0]+swy) gt !d.y_size
        doSub = (doSubx || doSuby)
        if (doSub) then begin
            nsubx = doSubx ? !d.x_size-px[0] : swx
            nsuby = doSuby ? !d.y_size-py[0] : swy
            sub_a = a[0:nsubx-1,0:nsuby-1]
        endif
    endif else begin		;Scale window
	if keyword_set(aspect) then begin
		if f ge 1.0 then swy = swy / f else swx = swx * f
		endif		;aspect
	tv,poly_2d(bytscl(a),$	;Have to resample image
		[[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
		keyword_set(interp),swx,swy), $
		px[0],py[0]
	endelse			;window_scale
  endelse			;scalable pixels

mx = !d.n_colors-1		;Brightest color
colors = [mx,mx,mx,0,0,0]	;color vectors

end