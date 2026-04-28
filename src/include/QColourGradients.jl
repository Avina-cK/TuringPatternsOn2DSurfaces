using Plots
gm7 = cgrad([colorant"rgb(7, 141, 112)"
                    , colorant"rgb(38, 206, 170)"
                    , colorant"rgb(152, 232, 193)"
                    , colorant"rgb(255,255,255)"
                    , colorant"rgb(123, 173, 226)"
                    , colorant"rgb(80, 73, 204)"
                    , colorant"rgb(61, 26, 120)"]
                );
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

nb4_bpyw_cat = cgrad([colorant"rgb(44, 44, 44)"
                , colorant"rgb(156, 89, 209)"
                
                , colorant"rgb(252, 244, 52)"
                , colorant"rgb(255, 255, 255)"
            ]
            ,4,categorical=true
           )
