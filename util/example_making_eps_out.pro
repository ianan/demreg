pro example_making_eps_out

  ; Example script just to show how to save nice .eps files out
  ;
  ;
  ; 18-Sep-2019 IGH   - Started
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ;

  ; Load in the batch file to set things up
  @post_outset
  
  ; as a single plot
  !p.multi=0
  ; But maybe want 2 rows and 3 columns
  ; !p.multi=[0,3,2]

  set_plot,'ps'
  device, /encapsulated, /color, /isolatin1, /inches, $
    bits=8, xsize=5, ysize=5,file='example_eps.eps'

  loadct,39,/silent
  
  plot_image,findgen(42,42),xtitle='x axis',ytitle='y axis',title='Title'

  device,/close
  set_plot, mydevice


  stop
end