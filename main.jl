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

ind = readtable("C:/PSID/csv/IND2011ER.csv")
tmp = readtable("C:/PSID/csv/FAM1968.csv")

yind = ind[[:ER30001, :ER30002, individual_relationship_to_head[1:length(years)], :ER30004, :ER30019]]
n = size(yind, 1)
yind = yind[yind[:ER30001].<=2930, :]
rename!(yind, :ER30001, :int_no_68)
rename!(yind, :ER30002, :person_no)
@printf "The original dataset has %d observations.\n" n
@printf "Dropping non-core individuals leaves %d observations\n" size(yind,1)

function create_panel()
  for yr = 1:length(years)
    @printf "\tWorking on data for year %d\n" years[yr]
    @time tmp = readtable("C:/PSID/csv/FAM"*string(years[yr])*".csv")
    @time rename!(tmp, family_1968_interview_number[yr], :int_no_68)
    #=@time yind = join(yind, tmp, on=:int_no_68)
    @time rename!(yind, individual_relationship_to_head[yr], convert(Symbol, "relhd_"*string(years[yr])))
    @time rename!(yind, family_head_annual_hours_prior_year[yr], convert(Symbol, "hd_hrs_"*string(years[yr])))
    @time rename!(yind, family_head_sex[yr], convert(Symbol, "hd_sex_"*string(years[yr])))
    @time rename!(yind, family_head_age[yr], convert(Symbol, "hd_age_"*string(years[yr])))
    @time rename!(yind, family_head_total_labor_income[yr], convert(Symbol, "hd_linc_"*string(years[yr])))=#
  end
end
