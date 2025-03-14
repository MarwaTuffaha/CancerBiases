function [comp_ref,comp_alt,comp_mer]=comp_3mer(ref,alt,mer)
    switch mer
        case 434
            comp_mer=121;
        case 334
            comp_mer=122;
        case 234
            comp_mer=123;
        case 134
            comp_mer=124;  
        case 433
            comp_mer=221;
        case 333
            comp_mer=222;
        case 233
            comp_mer=223;
        case 133
            comp_mer=224;
        case 432
            comp_mer=321;
        case 332
            comp_mer=322;  
        case 232
            comp_mer=323;
        case 132
            comp_mer=324;
        case 431
            comp_mer=421;
        case 331
            comp_mer=422;  
        case 231
            comp_mer=423;
        case 131
            comp_mer=424;            
        case 414
            comp_mer=141;
        case 314
            comp_mer=142;
        case 214
            comp_mer=143;
        case 114
            comp_mer=144;  
        case 413
            comp_mer=241;
        case 313
            comp_mer=242;
        case 213
            comp_mer=243;
        case 113
            comp_mer=244;
        case 412
            comp_mer=341;
        case 312
            comp_mer=342;  
        case 212
            comp_mer=343;
        case 112
            comp_mer=344;
        case 411
            comp_mer=441;
        case 311
            comp_mer=442;  
        case 211
            comp_mer=443;
        case 111
            comp_mer=444;
    end
    switch alt
        case 1
            comp_alt=4;
        case 2
            comp_alt=3;
        case 3
            comp_alt=2;
        case 4
            comp_alt=1;
    end    
    switch ref
        case 1
            comp_ref=4;
        case 2
            comp_ref=3;
        case 3
            comp_ref=2;
        case 4
            comp_ref=1;
    end     
end