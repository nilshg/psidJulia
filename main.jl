using DataFrames

include("varnames.jl")

(individual_age, individual_grades_completed,
 individual_id_interview_number, individual_relationship_to_head,
 individual_sample_weight, individual_sequence_number,
 family_1968_interview_number, family_head_age,
 family_head_annual_hours_prior_year, family_head_total_labor_income,
 family_head_sex, family_interview_number_by_year) = getvarnames()

years = [[1968:1997],[1999:2:2011]]
years = [1968:1993]


function create_panel()
  @time yind = readtable("C:/PSID/csv/IND2011ER.csv")
  yind = yind[[:ER30001, :ER30002, individual_relationship_to_head[1:length(years)],
               individual_age[1:length(years)]]]

  n = size(yind, 1)
  yind = yind[yind[:ER30001].<=2, :]
  rename!(yind, :ER30001, :int_no_68)
  rename!(yind, :ER30002, :person_no)
  for yr = 1:length(years)
    rename!(yind, individual_relationship_to_head[yr],
            convert(Symbol, "ind_relhd_"*string(years[yr])))
    rename!(yind, individual_age[yr],
            convert(Symbol, "ind_age_"*string(years[yr])))
  end
  @printf "The original dataset has %d observations.\n" n
  @printf "Dropping non-core individuals leaves %d observations\n" size(yind,1)
  yind

  for yr = 1:5#length(years)
    @printf "\tWorking on data for year %d\n" years[yr]
    # read in family data
    tmp = readtable("C:/PSID/csv/FAM"*string(years[yr])*".csv")
    # rename 1968 interview number
    rename!(tmp, family_1968_interview_number[yr], :int_no_68)
    # drop non-core households
    tmp = tmp[tmp[:int_no_68].<=2,:]
    rename!(tmp, family_interview_number_by_year[yr],
            convert(Symbol, "int_no_"*string(years[yr])))
    # only keep necessary columns
    tmp =
      tmp[[:int_no_68, :int_no_1969, family_head_age[yr], family_head_annual_hours_prior_year[yr],
             family_head_total_labor_income[yr], family_head_sex[yr]]]

    @printf "Before merge panel dimension: %d x %d\n" size(yind,1) size(yind,2)
    @printf "Dimension of remaining family panel: %d x %d\n" size(tmp,1) size(tmp,2)
    # merge
    yind = join(yind, tmp, on=:int_no_68, kind=:left)
    @printf "After merge panel dimension: %d x %d\n" size(yind,1) size(yind,2)
    # rename remaining columns
    rename!(yind, family_head_annual_hours_prior_year[yr],
              convert(Symbol, "hd_hrs_"*string(years[yr])))
    rename!(yind, family_head_sex[yr],
              convert(Symbol, "hd_sex_"*string(years[yr])))
    rename!(yind, family_head_age[yr],
              convert(Symbol, "hd_age_"*string(years[yr])))
    rename!(yind, family_head_total_labor_income[yr],
              convert(Symbol, "hd_linc_"*string(years[yr])))
  end

  return yind
end
