function getvarnames(symbols=true)

  individual_age =
    ["ER30004","ER30023","ER30046","ER30070","ER30094","ER30120","ER30141",
     "ER30163","ER30191","ER30220","ER30249","ER30286","ER30316","ER30346",
     "ER30376","ER30402","ER30432","ER30466","ER30501","ER30538","ER30573",
     "ER30609","ER30645","ER30692","ER30736","ER30809","ER33104","ER33204",
     "ER33304","ER33404","ER33504","ER33604","ER33704","ER33804","ER33904",
     "ER34004","ER34104"]

  individual_grades_completed =
    ["ER30010",  "NaN"  ,"ER30052","ER30076","ER30100","ER30126","ER30147",
     "ER30169","ER30197","ER30226","ER30255","ER30296","ER30326","ER30356",
     "ER30384","ER30413","ER30443","ER30478","ER30513","ER30549","ER30584",
     "ER30620","ER30657","ER30703","ER30748","ER30820","ER33115","ER33215",
     "ER33315","ER33415","ER33516","ER33616","ER33716","ER33817","ER33917",
     "ER34020","ER34119"]

  individual_id_interview_number =
    ["ER30001","ER30020","ER30043","ER30067","ER30091","ER30117","ER30138",
     "ER30160","ER30188","ER30217","ER30246","ER30283","ER30313","ER30343",
     "ER30373","ER30399","ER30429","ER30463","ER30498","ER30535","ER30570",
     "ER30606","ER30642","ER30689","ER30733","ER30806","ER33101","ER33201",
     "ER33301","ER33401","ER33501","ER33601","ER33701","ER33801","ER33901",
     "ER34001","ER34101"]

  individual_head_indicator =
    [[1 for i=1:15], [10 for i = 1:22]]

  individual_relationship_to_head =
    ["ER30003","ER30022","ER30045","ER30069","ER30093","ER30119","ER30140",
     "ER30162","ER30190","ER30219","ER30248","ER30285","ER30315","ER30345",
     "ER30375","ER30401","ER30431","ER30465","ER30500","ER30537","ER30572",
     "ER30608","ER30644","ER30691","ER30735","ER30808","ER33103","ER33203",
     "ER33303","ER33403","ER33503","ER33603","ER33703","ER33803","ER33903",
     "ER34003","ER34103"]

  individual_sample_weight =
    ["ER30019","ER30042","ER30066","ER30090","ER30116","ER30137","ER30159",
     "ER30187","ER30216","ER30245","ER30282","ER30312","ER30342","ER30372",
     "ER30398","ER30428","ER30462","ER30497","ER30534","ER30569","ER30605",
     "ER30641","ER30686","ER30730","ER30803","ER30864","ER33119","ER33275",
     "ER33318","ER33430","ER33546","ER33637","ER33740","ER33848","ER33950",
     "ER34045","ER34154"]

  individual_sequence_number =
    [  "NaN"  ,"ER30021","ER30044","ER30068","ER30092","ER30118","ER30139",
     "ER30161","ER30189","ER30218","ER30247","ER30284","ER30314","ER30344",
     "ER30374","ER30400","ER30430","ER30464","ER30499","ER30536","ER30571",
     "ER30607","ER30643","ER30690","ER30734","ER30807","ER33102","ER33202",
     "ER33302","ER33402","ER33502","ER33602","ER33702","ER33802","ER33902",
     "ER34002","ER34102"]

  family_1968_interview_number =
    ["V3",	"V534","V1230","V1932","V2533","V3085","V3497","V3909","V4423",
     "V5336","V5835","V6446","V7050","V7642","V8335","V8943","V10400","V11581",
     "V12988","V14090","V15105","V16605","V18021","V19321","V20621","V22400",
     "NaN","NaN","NaN","NaN","ER13019","ER17022","ER21009","ER25009","ER36009",
     "ER42009","ER47309"]

  family_head_age =
    ["V117","V1008","V1239","V1942","V2542","V3095","V3508","V3921","V4436",
     "V5350","V5850","V6462","V7067","V7658","V8352","V8961","V10419","V11606",
     "V13011","V14114","V15130","V16631","V18049","V19349","V20651","V22406"
    ,"ER2007","ER5006","ER7006","ER10009","ER13010","ER17013","ER21017",
     "ER25017","ER36017","ER42017","ER47317"]

  family_head_annual_hours_prior_year =
    [  "V47"  , "V465"  , "V1138" , "V1839" , "V2439" , "V3027"  ,"V3423" ,
     "V3823"  , "V4332" , "V5232" , "V5731" , "V6336" , "V6934" , "V7530" ,
     "V8228"  , "V8830" , "V10037", "V11146", "V12545", "V13745", "V14835",
     "V16335" , "V17744", "V19044", "V20344", "V21634", "ER4096",	"ER6936",
     "ER9187" ,"ER12174","ER16471","ER20399","ER24080","ER27886","ER40876",
     "ER46767","ER52175"]

  family_head_total_labor_income =
    ["V74","V514","V1196","V1897","V2498","V3051","V3463","V3863","V5031",
     "V5627","V6174","V6767","V7413","V8066","V8690","V9376","V11023","V12372",
     "V13624","V14671","V16145","V17534","V18878","V20178","V21484","V23323",
     "ER4140","ER6980","ER9231","ER12080","ER16463","ER20443","ER24116",
     "ER27931","ER40921","ER46829","ER52237"]

  family_head_sex =
    [   "V119",  "V1010",  "V1240",  "V1943",  "V2543",  "V3096",  "V3509",
       "V3922",  "V4437",  "V5351",  "V5851",  "V6463",  "V7068",  "V7659",
       "V8353",  "V8962", "V10420", "V11607", "V13012", "V14115", "V15131",
      "V16632", "V18050", "V19350", "V20652", "V22407",	"ER2008", "ER5007",
      "ER7007","ER10010","ER13011","ER17014","ER21018","ER25018","ER36018",
     "ER42018","ER47318"]

  family_interview_number_by_year =
    ["V3"  ,"V442","V1102","V1802","V2402","V3002","V3402","V3802","V4302",
     "V5202","V5702","V6302","V6902","V7502","V8202","V8802","V10002",
     "V11102","V12502","V13702","V14802","V16302","V17702","V19002","V20302",
     "V21602","ER2002","ER5002","ER7002","ER10002", "ER13002","ER17002",
     "ER21002","ER25002","ER36002","ER42002","ER47302"]

  family_education_head =
    ["V313",   "V794", "V1485", "V2197", "V2823", "V3241", "V3663", "V4198",
     "V5074",	"V5647", "V6194", "V6787", "V7433", "V8085", "V8709", "V9395",
     "V11042",	"V12400",	"V13640",	"V14687",	"V16161",	"V17545",	"V18898"]

  if symbols
    individual_age =
      [convert(Symbol, individual_age[i]) for i = 1:length(individual_age)]

    individual_grades_completed =
      [convert(Symbol, individual_grades_completed[i])
        for i = 1:length(individual_grades_completed)]

    individual_id_interview_number =
     [convert(Symbol, individual_id_interview_number[i])
       for i = 1:length(individual_id_interview_number)]

    individual_relationship_to_head =
      [convert(Symbol, individual_relationship_to_head[i])
        for i = 1:length(individual_relationship_to_head)]

    individual_sample_weight =
      [convert(Symbol, individual_sample_weight[i])
        for i = 1:length(individual_sample_weight)]

    individual_sequence_number =
      [convert(Symbol, individual_sequence_number[i])
        for i = 1:length(individual_sequence_number)]

    family_1968_interview_number =
      [convert(Symbol, family_1968_interview_number[i])
         for i = 1:length(family_1968_interview_number)]

    family_head_age =
      [convert(Symbol, family_head_age[i]) for i = 1:length(family_head_age)]

    family_head_annual_hours_prior_year =
      [convert(Symbol, family_head_annual_hours_prior_year[i])
         for i = 1:length(family_head_annual_hours_prior_year)]

    family_head_total_labor_income =
      [convert(Symbol, family_head_total_labor_income[i])
        for i = 1:length(family_head_total_labor_income)]

    family_head_sex =
      [convert(Symbol, family_head_sex[i]) for i = 1:length(family_head_sex)]

    family_interview_number_by_year =
      [convert(Symbol, family_interview_number_by_year[i])
         for i = 1:length(family_interview_number_by_year)]
  end

  return (individual_age, individual_grades_completed,
          individual_id_interview_number, individual_relationship_to_head,
          individual_sample_weight, individual_sequence_number,
          family_1968_interview_number, family_head_age,
          family_head_annual_hours_prior_year, family_head_total_labor_income,
          family_head_sex, family_interview_number_by_year)
end
