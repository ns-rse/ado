*! Program to generate indicators
*! Michael Hills 12/8/2002

program define genind
version 7.0
syntax [, EXPos(string) MOD(string) BASe(integer 1)]

if "`expos'" == "" {
    di in re "Exposure variable must be specified"
    exit
}
if "`mod'" == "" {
    di in re "Modifier variable must be specified"
    exit
}
qui inspect `expos'
if r(N_unique)>20 {
    di in re "Too many values for exposure"
    exit
}
qui inspect `mod'
if r(N_unique)>20 {
    di in re "Too many values for modifier"
    exit
}

local jlev 1
qui levels `mod' , local(modlev)
foreach j of local modlev {
    if `jlev' != `base' {
        gen byte I`mod'`jlev'=(`mod'==`j')
    } 
    local jlev=`jlev'+1
}

local ilev 1
qui levels `expos' , local(explev)
foreach i of local explev {
    local jlev 1
    foreach j of local modlev {
        if `ilev' != `base' {
            gen J`expos'`ilev'_`mod'`jlev'=(`expos'==`i')*(`mod'==`j')
        }
        local jlev=`jlev'+1
    }
    local ilev=`ilev'+1
}
end

    

