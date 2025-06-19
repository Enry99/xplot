xplot OUTCAR -r front --bonds none
xplot OUTCAR -r front --colorcode magmoms  --bonds none -o magmoms_cc.png
xplot OUTCAR -r front --bonds none --arrows magmoms --arrows-scale 1.5 -o magmoms_arrows.png