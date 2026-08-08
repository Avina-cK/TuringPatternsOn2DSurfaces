using Plots

gm7 = cgrad([colorant"rgb(7, 141, 112)"
                    , colorant"rgb(38, 206, 170)"
                    , colorant"rgb(152, 232, 193)"
                    , colorant"rgb(255,255,255)"
                    , colorant"rgb(123, 173, 226)"
                    , colorant"rgb(80, 73, 204)"
                    , colorant"rgb(61, 26, 120)"]
                );

les7 = cgrad([colorant"rgb(213, 45, 0)"
                    , colorant"rgb(239, 118, 39)"
                    , colorant"rgb(255, 154, 86)"
                    , colorant"rgb(255,255,255)"
                    , colorant"rgb(209, 98, 164)"
                    , colorant"rgb(181, 86, 144)"
                    , colorant"rgb(163, 2, 98)"]
                )

nb4_bpyw = cgrad([colorant"rgb(44, 44, 44)"
                , colorant"rgb(156, 89, 209)"
                , colorant"rgb(252, 244, 52)"
                , colorant"rgb(255, 255, 255)"
            ]
           )

rainbow6 = cgrad([
                colorant"rgb(228,3,3)",
                colorant"rgb(255,140,0)",
                colorant"rgb(255,237,0)",
                colorant"rgb(0,128,38)",
                colorant"rgb(0,76,255)",
                colorant"rgb(115,41,130)"
                ])

tragen3 = cgrad([
                colorant"rgb(91,206,250)",
                colorant"rgb(255,255,255)",
                colorant"rgb(245,169,184)"
                ])

bis3 = cgrad([
            colorant"rgb(214,2,112)",
            colorant"rgb(155,79,150)",
            colorant"rgb(0,56,168)"
            ])

ace4 = cgrad([
            colorant"rgb(0,0,0)",
            colorant"rgb(163,163,163)",
            colorant"rgb(255,255,255)",
            colorant"rgb(128,0,128)"
            ])

pan3 = cgrad([
            colorant"rgb(255,33,140)",
            colorant"rgb(255,216,0)",
            colorant"rgb(33,177,255)"
            ])

aro5 = cgrad([
            colorant"rgb(61,165,66)",
            colorant"rgb(167,211,121)",
            colorant"rgb(255,255,255)",
            colorant"rgb(169,169,169)",
            colorant"rgb(0,0,0)"
])

aroace5 = cgrad([
            colorant"rgb(226, 140, 0)",
            colorant"rgb(236, 205, 0)",
            colorant"rgb(255,255,255)",
            colorant"rgb(98, 174, 220)",
            colorant"rgb(32, 56, 86)"
])

genque3 = cgrad([
            colorant"rgb(181,126,220)",
            colorant"rgb(255,255,255)",
            colorant"rgb(74,129,35)"
            ])

gm7_r = cgrad([colorant"rgb(61, 26, 120)"
                    , colorant"rgb(80, 73, 204)"
                    , colorant"rgb(123, 173, 226)"
                    , colorant"rgb(255,255,255)"
                    , colorant"rgb(152, 232, 193)"
                    , colorant"rgb(38, 206, 170)"
                    ,colorant"rgb(7, 141, 112)"]
                );
                
gm7_cat_r = cgrad([colorant"rgb(61, 26, 120)"
                    , colorant"rgb(80, 73, 204)"
                    , colorant"rgb(123, 173, 226)"
                    , colorant"rgb(255,255,255)"
                    , colorant"rgb(152, 232, 193)"
                    , colorant"rgb(38, 206, 170)"
                    ,colorant"rgb(7, 141, 112)"]
                , 7
                , categorical=true
                );
gm7_cat = cgrad([colorant"rgb(7, 141, 112)"
                    , colorant"rgb(38, 206, 170)"
                    , colorant"rgb(152, 232, 193)"                    
                    , colorant"rgb(255,255,255)"
                    , colorant"rgb(123, 173, 226)"
                    , colorant"rgb(80, 73, 204)"
                    , colorant"rgb(61, 26, 120)"]
                , 7
                , categorical=true
                );


les9_cat = cgrad([colorant"rgb(213, 45, 0)"
                    , colorant"rgb(239, 118, 39)"
                    , colorant"rgb(255, 154, 86)"
                    , colorant"rgb(255,255,255)"
                    , colorant"rgb(209, 98, 164)"
                    , colorant"rgb(181, 86, 144)"
                    , colorant"rgb(163, 2, 98)"]
                , 9, categorical=true)

                
nb4_bpyw_cat = cgrad([colorant"rgb(44, 44, 44)"
                , colorant"rgb(156, 89, 209)"
                
                , colorant"rgb(252, 244, 52)"
                , colorant"rgb(255, 255, 255)"
            ]
            ,4,categorical=true
           )
