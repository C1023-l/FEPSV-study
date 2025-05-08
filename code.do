cd "*******"
import delimited "*******", clear

* 使用logit模型估计倾向得分

ssc install psmatch2, replace
logit treated IC PD RL BUS GSA PGDP
predict ps_score, pr
psmatch2 treated (IC PD RL BUS GSA PGDP), out (ps_score) neighbor(3) caliper(0.05) logit

* 将匹配结果保存
gen matched = _weight > 0
pstest IC PD RL BUS GSA PGDP, graph
pstest IC PD RL BUS GSA PGDP, detail

clear
*（1）生成交互项---基准回归数据集构造
cd "*******"
use "*******"
joinby "*******"
gen time=0
replace time=1 if date>318//date=2023年11月15日，为2023年的第318天。设置对照组，并标记为1
gen did= treated * time 

*（2）基准回归
use "*******"
gen log_AQI=ln(aqi) //因变量取对数
duplicates drop city_code date,force 
xtset city_code date 
xtreg log_AQI did IC PD RL BUS GSA PGDP 降水 风速 温度 i.date,fe //面板回归并控制个体固定效应、时间效应
xtreg log_AQI did i.date, fe //不考虑控制变量

*（3）平行趋势检验
gen xianhou=date-startdate if treated1111>0 
tab xianhou 
tab xianhou, gen(xh)
forvalues i=1/250{
replace xh`i'=0 if xh`i'==.
} //将所有空值设置为0
drop xh1
xtreg log_AQI did IC PD RL BUS GSA PGDP 降水 风速 温度 i.date,fe 
coefplot, baselevels keep(xh*) vertical yline(0) ytitle("政策效应") xtitle("政策前天数")addplot(line @b @at) ciopts(recast(rcap)) scheme(s1mono) levels(95) coeflabels(xh2 = "-250" xh3 = "-200" xh4 = "-150" xh5 = "-100" xh6 = "-50") xline(6,lp(shortdash))//绘图命令

*（4）安慰剂检验
clear
set matsize 1100
mat b=J(1000,1,.)
mat se=J(1000,1,.)
mat p=J(1000,1,.)
cd *********** 

forvalues i= 1/1000{
 use ******,clear
 duplicates drop city_code date,force
 xtset city_code date
 keep if date==318
 sample 15,count 
 keep city_code
 save matchid.dta,replace
 merge 1:m city_code using *******
 gen treat=(_merge==3)
 save matchid`i'.dta,replace
 
 use *******
 bsample 1, strata(city_code)
 keep date
 save matchdate.dta, replace
 mkmat date,matrix(sampledate)
 
 use matchid`i'.dta,replace
 duplicates drop city_code date,force
 xtset city_code date
 gen time=0
 foreach j of numlist 1/285 {
 replace time = 1 if (city_code == `j'& date>=sampledate[`j',1])
 }
 gen did1=time*treat
 qui xtreg log_AQI did IC PD RL BUS GSA PGDP 降水 风速 温度 i.date,fe
 mat b[`i',1]=_b[did1]
 mat se[`i',1]=_se[did1]
 scalar df_r=e(N)-e(df_m)-1
 mat p[`i',1]=2*ttail(df_r,abs(_b[did1]/_se[did1]))
 }

svmat b, names(coef)
svmat se, names(se)
svmat p, names(pvalue)
drop if pvalue1 ==.
label var pvalue1 p值
label var coef1 估计系数

preserve
    kdensity coef1, generate(density_x density_y) nograph
    save "temp_density.dta", replace
restore

* 载入核密度数据
merge 1:1 _n using "temp_density.dta", nogen keepusing(density_x density_y)

twoway ///
    (scatter pvalue1 coef1, msymbol(O) mcolor(red)) ///
    (line density_y density_x, lcolor(brown) lwidth(medthick) yaxis(2)), ///
    xlabel(-0.5(0.1)0.5, grid) ///
    ytitle("p值", axis(1)) ///
    ytitle("核密度", axis(2)) ///
    xtitle("估计系数") ///
    yline(0.1, lp(shortdash)) ///
	xline(0.4228,lp(shortdash)) ///
    legend(off) ///
    title("安慰剂检验")

	
*（4）稳健性检验之反事实检验
*生成交互项---基准回归数据集构造
cd "*******" //打开文件夹

use "*******"
gen advance_100=data-100
gen advance_200=date-200
gen forward_100=date+100
gen forward_200=date+200
save "*******"

/*cd ****** //打开文件夹
use ****** //打开所用数据
joinby *********
replace time=1 if date>318//date=2023年11月15日，为2023年的第318天。设置对照组，并标记为o
gen did= treated* time //生成DID交互项

use ***********
gen log_AQI=ln(aqi) //因变量取对数
duplicates drop city_code date,force //删除重复值
xtset city_code date //设置面板数据
xtreg log_AQI did IC PD RL BUS GSA PGDP 降水 风速 温度 i.date,fe //面板回归并控制个体固定效应、时间效应
xtreg log_AQI did i.date, fe //不考虑控制变量*/

*反事实检验之提前100天----交互项构造
replace advance_100=2022 if treated1111==0 //提前100天——将对照组开始日期设置为“2022”
gen t_advance_100=0
replace t_advance_100=1 if date>= advance_100
gen didadvance_100= treated1111* t_advance_100 

*反事实检验之提前200天----交互项构造
replace advance_200=2022 if treated1111==0 //提前100天——将对照组开始日期设置为“2022”
gen t_advance_200=0
replace t_advance_200=1 if date>= advance_200 //生成政策变量
gen didadvance_200= treated1111* t_advance_200 //生成DID交互项

*反事实检验之推后100天----交互项构造
replace forward_100=2022 if treated1111==0 //提前100天——将对照组开始日期设置为“2022”
gen t_forward_100=0
replace t_forward_100=1 if date>= forward_100 //生成政策变量
gen didforward_100= treated1111* t_forward_100 //生成DID交互项

*反事实检验之推后200天----交互项构造
replace forward_200=2022 if treated1111==0 //提前100天——将对照组开始日期设置为“2022”
gen t_forward_200=0
replace t_forward_200=1 if date>= forward_200 //生成政策变量
gen didforward_200= treated1111* t_forward_200 //生成DID交互项

drop advance_100 advance_200 forward_100 forward_200
drop t_advance_100 t_advance_200 t_forward_100 t_forward_200

save "*********" //保存为最终的数据
*面板数据设定

gen log_commuting=ln(commuting) //因变量取对数
duplicates drop city_code date,force //删除重复值
xtset city_code date //设置面板数据

*回归
*提前100天
xtreg log_AQI didadvance_100 IC PD RL BUS GSA PGDP 降水 风速 温度 i.date,fe //面板回归并控制个体固定效应、时间效应
*提前200天
xtreg log_AQI didadvance_200 IC PD RL BUS GSA PGDP 降水 风速 温度 i.date,fe //面板回归并控制个体固定效应、时间效应
*推后100天
xtreg log_AQI didforward_200 IC PD RL BUS GSA PGDP 降水 风速 温度 i.date,fe //面板回归并控制个体固定效应、时间效应
*推后200天
xtreg log_AQI didforward_100 IC PD RL BUS GSA PGDP 降水 风速 温度 i.date,fe //面板回归并控制个体固定效应、时间效应