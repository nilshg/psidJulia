using DataFrames

function prepCovmat()

  # Some parameters
  data = readtable("C:/Users/tew207/Desktop/ready_newdata.csv")
  tinit = 67; tlast = 96
  ageinit = 20; agelast = 64
  agecell = 4
  minyrs = 20
  nlag = 29
  agelb = 19; ageub = agelb + agecell
  agemidpt = (agelb+ageub)/2
  agemax = (agelast + agelast - agecell)/2 - (ageub+agelb)/2
  oldcoh = agelast - minyrs - ageinit
  newcoh = tlast - (minyrs-1) - (tinit+1)
  maxcoh = oldcoh + newcoh
  for i = 94:97
    rename!(data, convert(Symbol, "upedu"*string(i)*"h"),
                  convert(Symbol, "grade"*string(i)))
  end

  #####################
  ###   EDUCATION  ####
  #####################

  data[:grade72] = data[:edcn72]
  data[array(data[:grade72].>25, 0), :] = NA # take out values of grade above 25
  data[:grade75] = data[:edcn75]
  data[array(data[:grade75].>25, 0), :] = NA

  @where(data, array((:seqno68 .== 1) & (:seqno69 .== 1),false))[:grade69]
    = @where(data, array((:seqno68 .== 1) & (:seqno69 .== 1),false))data[:grade68]

  @where(data, array((:grade69.==0) & (:seqno69.==1) & (:seqno72.==1), false))[:grade69]
    = @where(data, array((:grade69.==0) & (:seqno69.==1) & (:seqno72.==1), false))[:grade72]

  @where(data , array((:seqno68.==1) & (:seqno70.==1), false))[:grade70] =
    @where(data , array((:seqno68.==1) & (:seqno70.==1), false))[:grade68]

  @where(data, array((:grade69).==0 & (:seqno70.==1) & (:seqno72.==1), false))[:grade70] =
    @where(data, array((:grade69).==0 & (:seqno70.==1) & (:seqno72.==1), false))[:grade72]

  @where(data, array((:seqno70 .== 1) & (:seqno71 .== 1),false))[:grade71] =
    @where(data, array((:seqno70 .== 1) & (:seqno71 .== 1),false))[:grade70]

  [:grade71] = [:grade70]  (:seqno70.==1) & (:seqno71.==1)
  [:grade71] = [:grade72]  (:grade71).==0 & (:seqno71.==1) & (:seqno72.==1)

  [:grade73] = [:grade72]  (:seqno72.==1) & (:seqno73.==1)
  [:grade73] = [:grade75]  (:grade73).==0 & (:seqno73.==1) & (:seqno75.==1)

  [:grade74] = [:grade72]  (:seqno72.==1) & (:seqno74.==1)
  [:grade74] = [:grade75]  (:grade74).==0 & (:seqno74.==1) & (:seqno75.==1)

  for i = (tinit+1):(tlast+1)  # take out grades above 30 (uncertain!)
    data[data[convert(Symbol, "grade"*string(i))].>30, :] = NA
  end

  #################
  ###   WAGES   ###
  #################

  awg = {"67"=>2.85, "68"=>3.02, "69"=>3.22, "70"=>3.40, "71"=>3.63, "72"=>3.90,
         "73"=>4.14, "74"=>4.43, "75"=>4.73, "76"=>5.06, "77"=>5.44, "78"=>5.87,
         "79"=>6.33, "80"=>6.84, "81"=>7.43, "82"=>7.86, "83"=>8.19, "84"=>8.48,
         "85"=>8.73, "86"=>8.92, "87"=>9.13, "88"=>9.43, "89"=>9.80, "90"=>10.19,
         "91"=>10.50, "92"=>10.76, "93"=>11.03, "94"=>11.32, "95"=>11.64,
         "96"=>12.03}

  for i = tinit:tlast
    data[convert(Symbol, "rawg"*string(i))] =
      awg[string(i)]./data[convert(Symbol, "prc"*string(i))]

    @where(data, array((convert(Symbol, "rhdwg"*string(i)).<=2*convert(Symbol, "rawg"*string(i))/awg["96"]) |
                       (convert(Symbol, "rhdwg"*string(i)).>400*convert(Symbol, "rawg"*string(i))/awg["96"]) |
                       (convert(Symbol, "hwkhrs"*string(i)).> 5096) | (convert(Symbol, "hwkhrs"*string(i)).<520) ,false))
  end
end
