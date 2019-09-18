; Batch script to make sure default IDL procedure plots are nice
;
; 18-Sep-2019 IGH Started
;
;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

clearplot
!y.range=0
!x.range=0
!p.title=''
!x.title=''
!y.title=''

; For x and y range to exactly values of data
!y.style=17
!x.style=17

; if on windows then this should be set_plot,'win'
set_plot,'x'
device,retain=2, decomposed=0
mydevice = !d.name

;; Make IDL use device/hardware fonts
!p.font = 0

;; Other things to make graphs nicer
!p.color = 255  ;255 for white line
!p.background = 0   ;0 for balck background
!p.thick = 1
!p.charthick = 1
!p.charsize = 1.
!p.symsize=1

;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
