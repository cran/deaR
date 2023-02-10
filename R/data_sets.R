#' Data: Charnes, Cooper and Rhodes (1981).
#'
#' Data from Project Follow Through (PTF) in public school education. There are
#' 49 DMUs (school sites) in PFT and  21 DMUs in Non-Follow Through (NFT). Authors
#' consider 3 outputs (Y) and 5 inputs (X).
#' 
#' @usage data("PFT1981")
#' 
#' @format Data frame with 70 rows and 10 columns. Definition of inputs (X) and
#' outputs (Y):
#' \describe{
#'   \item{Y1 = Reading}{Total Reading Scores (as measured by the Metropolitan
#'   Achievement Test).}
#'   \item{Y2 = Math}{Total Math Scores (total mathematics score by the
#'   Metropolitan Achievement Test.}
#'   \item{Y3 = Coopersmith}{Total Coopersmith Scores (Coopersmith self-esteem
#'   inventory, intended as a measure of self-esteem).}
#'   \item{X1 = Education}{Education level of mother (as measured in terms of
#'   percentage of high school graduates among female parents).}
#'   \item{X2 = Occupation}{Occupation Index (highest occupation of a family
#'   member according to a pre-arranged rating scale).}
#'   \item{X3 = Parental}{Parental Visit Index (representing the number of
#'   visits to the school site).}
#'   \item{X4 = Counseling}{Counseling Index (parent counselling index
#'   calculated from data on time spent with child on school-related topics such
#'   as reading together, etc.).}
#'   \item{X5 = Teachers}{Number of Teachers (number of teachers at a given site).}
#'   \item{Program}{PFT or NFT.}
#' }
#' 
#' @source Charnes, A.; Cooper, W.W.; Rhodes, E. (1981). "Evaluating Program and
#' Managerial Efficiency: An Application of Data Envelopment Analysis to Program
#' Follow Through", Management Science, 27(6), 668-697.
#' \doi{10.1287/mnsc.27.6.668}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example 1. Replication of results in Charnes, Cooper and Rhodes (1981)
#' data("PFT1981")
#' # selecting DMUs in Project Follow Through (PFT)
#' PFT <- PFT1981[1:49, ]
#' PFT <- make_deadata(PFT,
#'                     dmus = 1,
#'                     inputs = 2:6,
#'                     outputs = 7:9 )
#' eval_pft <- model_basic(PFT,
#'                         orientation = "io",
#'                         rts = "crs")
#' eff_pft <- efficiencies(eval_pft)
#'
#' # Example 2. Replication of results in Charnes, Cooper and Rhodes (1981)
#' data("PFT1981")
#' # selecting DMUs in Non-Follow Through (NFT)
#' NFT <- PFT1981[50:70,]
#' NFT <- make_deadata(NFT,
#'                     dmus = 1,
#'                     inputs = 2:6,
#'                     outputs = 7:9 )
#' eval_nft <- model_basic(NFT,
#'                         orientation = "io",
#'                         rts = "crs")
#' eff_nft <- efficiencies(eval_nft)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_basic}}

"PFT1981"


#' Data: Tone (2002).
#'
#' This dataset consists of six power plants with 4 inputs (X) and 2 outputs (Y).
#' 
#' @usage data("Power_plants")
#' 
#' @format Data frame with 15 rows and 7 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1}{Manpower requiered}
#'   \item{x2}{Construction costs in millions of dollars}
#'   \item{x3}{Annual maintenance costs in millions of dollars}
#'   \item{x4}{Number of villages to be evacuated }
#'   \item{y1}{Power generated in megawatts}
#'   \item{y2}{Safety level}
#' }
#'
#' @source Andersen, P.; Petersen, N.C. (1993). "A procedure for ranking efficient
#' units in data envelopment analysis", Management Science, 39, 1261-1264.
#'
#' Doyle, J. and Green R. (1993). "Data envelopment analysis and multiple criteria
#' decision making", Omega, 21 (6), 713-715.  \doi{10.1016/0305-0483(93)90013-B}
#'
#' Tone, K. (2002). "A slacks-based measure of super-efficiency in data envelopment
#' analysis", European Journal of Operational Research, 143, 32-41.
#' \doi{10.1016/S0377-2217(01)00324-1}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example 1. Radial super-efficiency model.
#' # Replication of results in Tone (2002)
#' data("Power_plants")
#' data_example <- make_deadata(Power_plants,
#'                              ni = 4,
#'                              no = 2)
#' result <- model_supereff(data_example,
#'                          orientation = "io",
#'                          rts = "crs")
#' eff <- efficiencies(result)
#' eff
#'
#' # Example 2. SBM super-efficiency model.
#' data("Power_plants")
#' data_example <- make_deadata(Power_plants,
#'                              ni = 4,
#'                              no = 2)
#' result2 <- model_sbmsupereff(data_example,
#'                              orientation = "io",
#'                              rts = "crs")
#' efficiencies(result2)
#' slacks(result2)$input
#' references(result2)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_supereff}},
#' \code{\link{model_sbmsupereff}}

"Power_plants"


#' Data: Zhu (2014).
#'
#' This dataset consists of 15 firms from the Fortune 500 list 1995
#' (\url{https://fortune.com/ranking/fortune500/}) with 3 inputs and 2 outputs.
#' 
#' @usage data("Fortune500")
#' 
#' @format Data frame with 15 rows and 6 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1 = Assets}{Assets (millions of dollars)}
#'   \item{x2 = Equity}{Equity (millions of dollars)}
#'   \item{x3 = Employees}{Number of employees}
#'   \item{y1 = Revenue}{Revenue (millions of dollars)}
#'   \item{y2 = Profit}{Profit (millions of dollars)}
#' }
#' 
#' @source Zhu, J. (2014). Quantitative Models for Performance Evaluation and
#' Benchmarking. Data Envelopment Analysis with Spreadsheets. 3rd Edition Springer,
#' New York. \doi{10.1007/978-3-319-06647-9}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#' 
#' \strong{Vicente Bolos} (\email{vicente.bolos@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' data("Fortune500")
#' data_Fortune <- make_deadata(datadea = Fortune500,
#'                              dmus = 1,
#'                              inputs = 2:4,
#'                              outputs = 5:6)
#' result <- model_multiplier(data_Fortune,
#'                            epsilon = 1e-6,
#'                            orientation = "io",
#'                            rts = "crs")
#' # results for General Motors and Ford Motor are not shown
#' # by deaR because the solution is infeasible
#' efficiencies(result)
#' multipliers(result)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_multiplier}}

"Fortune500"


#' Data: Wu, Tsai and Zhou (2011).
#'
#' This dataset consists of 23 four- and five-plum ITHs in Taipei in 2006. Authors
#' consider 4 inputs and 3 outputs.
#' 
#' @usage data("Hotels")
#' 
#' @format Data frame with 23 rows and 8 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1 = Employees}{Total number of employees)}
#'   \item{x2 = Guest_rooms}{Total number of guest rooms)}
#'   \item{x3 = Area_F&B}{Total area of F&B departments (in 36 square-feet)}
#'   \item{x4 = Operating_cost}{Total operating cost (in NT$)}
#'   \item{y1 = Room_revenue}{Room revenues (in NT$)}
#'   \item{y2 = F&B_revenue}{F&B revenues (in NT$)}
#'   \item{y3 = Other_revenue}{Other revenues (in NT$)}
#' }
#' 
#' @source Wu, J.; Tsai, H. and Zhou, Z. (2011). "Improving efficiency in International
#' tourist hotels in Taipei using a non-radial DEA mode", Internationl Journal of
#' Contemporary Hospitality Management, 23(1), 66-83. \doi{10.1108/09596111111101670}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. Replication of results in Wu,Tsai and Zhou (2011)
#' data("Hotels")
#' data_hotels <- make_deadata(Hotels,
#'                             dmus = 1,
#'                             inputs = 2:5,
#'                             outputs = 6:8)
#' result <- model_nonradial(data_hotels,
#'                           orientation = "oo",
#'                           rts = "vrs")
#' efficiencies(result)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_nonradial}}

"Hotels"


#' Data: Tomkins and Green (1988).
#'
#' Data from 20 University accounting departments in the UK.
#' 
#' @usage data("Departments")
#' 
#' @format Data frame with 20 rows and 11 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1 = Staff}{Average Full Time Academic Staff 82/3-84/5)}
#'   \item{x2 = Salaries}{1984-5 Salaries Academics and Related (in pounds))}
#'   \item{x3 = Other_Exp}{1984-5 Other Expenses (in pounds)}
#'   \item{y1 = Undergrad}{Average Number Undergraduates 82/3-84/5}
#'   \item{y2 = Research_post}{Research Postgraduates }
#'   \item{y3 = Taught_post}{Taught Postgraduates}
#'   \item{y4 = Res_co_income}{Research council income (in pounds)}
#'   \item{y5 = Other_res_income}{Other research income (in pounds)}
#'   \item{y6 = Other_income}{Other income (in pounds)}
#'   \item{y7 = Publications}{Number of publications}
#' }
#' 
#' @source Tomkins, C.; Green, R. (1988). "An Experiment in the Use of Data Envelopment
#' Analysis for Evaluating the Efficiency of UK University Departments of Accounting",
#' Financial Accountability and Management, 4(2), 147-164.
#' \doi{10.1111/j.1468-0408.1988.tb00296.x}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example.
#' # Replication of results DEA1 in Tomkins and Green (1988)
#' data("Departments")
#' # Calculate Total income
#' Departments$Total_income <- Departments[, 5] + Departments[, 6] + Departments[, 7]
#' data_example <- make_deadata(Departments,
#'                              inputs = 9,
#'                              outputs = c(2, 3, 4, 12))
#' result <- model_basic(data_example,
#'                       orientation = "io",
#'                       rts = "crs")
#' efficiencies(result) # Table 3 (p.156)
#' references(result) # Table 3 (p.157)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_basic}}

"Departments"


#' Data: Sanei and Mamizadeh Chatghayeb (2013).
#'
#' Data of 17 supply chain (buyer-supplier relationship in manufacturing).
#' 
#' @usage data("Supply_Chain")
#' 
#' @format Data frame with 17 rows and 8 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{X1 to X3}{Inputs of buyers}
#'   \item{I1 to I2}{Outputs of buyers, Inputs of suppliers}
#'   \item{Y1 to Y2}{Outputs of suppliers}
#' }
#' 
#' @source Sanei, M.; Mamizadeh Chatghayeb, S. (2013). “Free Disposal Hull Models
#' in Supply Chain Management”, International Journal of Mathematical Modelling
#' and Computations, 3(3), 125-129.
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. FDH input-oriented.
#' # Replication of results in Sanei and Mamizadeh Chatghayeb (2013)
#' data("Supply_Chain")
#' data_fdh1 <- make_deadata(Supply_Chain,
#'                           dmus = 1,
#'                           inputs = 2:4,
#'                           outputs = 5:6)
#' # by default orientation = "io"
#' result <- model_fdh(data_fdh1)
#' efficiencies(result)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_fdh}}

"Supply_Chain"


#' Data: Tone (2001).
#'
#' Data of 5 DMUs producing 2 outputs by using 2 inputs
#' 
#' @usage data("Tone2001")
#' 
#' @format Data frame with 5 rows and 5 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1}{Input1}
#'   \item{x2}{Input2}
#'   \item{y1}{Output1}
#'   \item{y2}{Output2}
#' }
#' 
#' @source Tone, K. (2001). "A Slacks-Based Measure of Efficiency in Data Envelopment
#' Analysis", European Journal of Operational Research, 130, 498-509.
#' \doi{10.1016/S0377-2217(99)00407-5}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. Replication of results in Tone (2001, p. 505)
#' data("Tone2001")
#' data_example <- make_deadata(Tone2001,
#'                              ni = 2,
#'                              no = 2)
#' result <- model_sbmeff(data_example,
#'                        orientation = "no",
#'                        rts = "crs")
#' efficiencies(result)
#' slacks(result)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_sbmeff}}

"Tone2001"


#' Data: Wang and Lan (2011).
#'
#' Data of the industrial economy of China in 2005-2009 (data in wide format).
#' 
#' @usage data("Economy")
#' 
#' @format Data frame with 31 rows and 16 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1 = Capital}{Total assets (in 100 million RMB)}
#'   \item{x2 = Labor}{Annual average employed persons (in 10000 persons)}
#'   \item{y1 = GIOV}{Gross industrial output value (in 100 million RMB)}
#' }
#' 
#' @source Wang, Y.; Lan, Y. (2011). "Measuring Malmquist Productiviy Index: A New
#' Approach Based on Double Frontiers Data Envelopment Analysis". Mathematical and
#' Computer Modelling, 54, 2760-2771. \doi{10.1016/j.mcm.2011.06.064}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example . Data in wide format.
#' # Replication of results in Wang and Lan (2011, p. 2768)
#' data("Economy")
#' data_example <- make_malmquist(Economy,
#'                                nper = 5,
#'                                arrangement = "horizontal",
#'                                ni = 2,
#'                                no = 1)
#' result <- malmquist_index(data_example)
#'
#' @seealso \code{\link{make_malmquist}}, \code{\link{malmquist_index}}

"Economy"


#' Data: Wang and Lan (2011).
#'
#' Data of the industrial economy of China in 2005-2009 (data in long format).
#' 
#' @usage data("EconomyLong")
#' 
#' @format Data frame with 155 rows and 5 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1 = Capital}{Total assets (in 100 million RMB)}
#'   \item{x2 = Labor}{Annual average employed persons (in 10000 persons)}
#'   \item{y1 = GIOV}{Gross industrial output value (in 100 million RMB)}
#' }
#' 
#' @source Wang, Y.; Lan, Y. (2011). "Measuring Malmquist Productiviy Index: A
#' New Approach Based on Double Frontiers Data Envelopment Analysis". Mathematical
#' and Computer Modelling, 54, 2760-2771. \doi{10.1016/j.mcm.2011.06.064}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. Data in long format.
#' # Replication of results in Wang and Lan (2011, p. 2768)
#' data("EconomyLong")
#' data_example <- make_malmquist(EconomyLong,
#'                                percol = 2,
#'                                arrangement = "vertical",
#'                                ni = 2,
#'                                no = 1)
#' result <- malmquist_index(data_example)
#'
#' @seealso \code{\link{make_malmquist}}, \code{\link{malmquist_index}}

"EconomyLong"


#' Data: Golany and Roll (1989).
#'
#' Data of 13 DMUs using 3 inputs to produce 2 outputs.
#' 
#' @usage data("Golany_Roll_1989")
#' 
#' @format Data frame with 13 rows and 6 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1}{Input 1}
#'   \item{x2}{Input 2}
#'   \item{x3}{Input 3}
#'   \item{y1}{Output 1}
#'   \item{y1}{Output 2}
#' }
#' 
#' @source Golany, B.; Roll, Y. (1989). "An Application Procedure for DEA". Omega,
#' International Journal of Management Science, 17(3), 237-250.
#' \doi{10.1016/0305-0483(89)90029-7}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example.
#' data("Golany_Roll_1989")
#' data_example <- make_deadata(datadea = Golany_Roll_1989,
#'                              dmus = 1,
#'                              inputs = 2:4,
#'                              outputs = 5:6)
#' result <- cross_efficiency(data_example,
#'                            orientation = "io",
#'                            selfapp = TRUE)
#' result$Arbitrary$cross_eff
#' result$Arbitrary$e
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_multiplier}}, \code{\link{cross_efficiency}}

"Golany_Roll_1989"


#' Data: Doyle and Green (1994).
#'
#' Data adapted from Tomkins and Green (1988).  13 DMUs using 3 inputs to produce 2 outputs.
#' 
#' @usage data("Doyle_Green_1994")
#' 
#' @format Data frame with 13 rows and 6 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{y1 = Undergraduate}{Number of undergraduates}
#'   \item{y2 = Postgraduates}{Number of postgraduates (taught and research)}
#'   \item{y3 = Research_income}{Research and other income}
#'   \item{y4 = Publications}{Number of publications}
#'   \item{x1 = Salaries}{Salaries of academic and related staff}
#'   \item{x2 = Other_expenses}{Other expenses}
#' }
#' 
#' @source Doyle, J.; Green, R. (1994). “Efficiency and cross efficiency in DEA:
#' derivations, meanings and the uses”,  Journal of Operational Research Society,
#' 45(5), 567–578. \doi{10.2307/2584392}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example.
#' data("Doyle_Green_1994")
#' data_example <- make_deadata(datadea = Doyle_Green_1994,
#'                             dmus = 1,
#'                             inputs = 6:7,
#'                             outputs = 2:5)
#' result <- cross_efficiency(data_example,
#'                            orientation = "io",
#'                            selfapp = TRUE)
#' result$Arbitrary$cross_eff
#' result$Arbitrary$e
#' # Aggressive using method II
#' result$M2_agg$cross_eff
#' # Aggressive using method III
#' result$M3_agg$cross_eff
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_multiplier}},
#' \code{\link{cross_efficiency}}

"Doyle_Green_1994"


#' Data: Cooper, Seiford and Tone (2007).
#'
#' Data for 23 public libraries of the Tokyo Metropolitan Area in 1986.
#' 
#' @usage data("Libraries")
#' 
#' @format Data frame with 23 rows and 7 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1 = AREA}{Floor area (unit=1000 m2)}
#'   \item{x2 = BOOKS}{Number of books (unit=1000)}
#'   \item{x3 = STAFF}{Staff}
#'   \item{x4 = POPULATION}{Population (unit=1000)}
#'   \item{y1 = REGISTERED}{Registered residents (unit=1000)}
#'   \item{y2 = BORROWED}{Borrowed books (unit=1000)}
#' }
#' 
#' @source Cooper, W.W.; Seiford, L.M. and Tone, K. (2007). Data Envelopment Analysis.
#' A Comprehensive Text with Models, Applications, References and DEA-Solver Software.
#' Springer.
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example 1. Non-controllable input (POPULATION).
#' # Replication of results in Cooper, Seiford and Tone (2007, p.221)
#' data(Libraries)
#' # POPULATION (non-controllable input) is the 4th input.
#' data_example <- make_deadata(Libraries,
#'                              dmus = 1,
#'                              inputs = 2:5,
#'                              nc_inputs = 4,
#'                              outputs = 6:7)
#' result <- model_basic(data_example,
#'                       orientation = "io",
#'                       rts = "crs")
#' efficiencies(result)
#' targets(result)
#'
#' # Example 2. Non-discretionary input (POPULATION).
#' data(Libraries)
#' # POPULATION (non-controllable input) is the 4th input.
#' data_example2 <- make_deadata(Libraries,
#'                               dmus=1,
#'                               inputs=2:5,
#'                               nd_inputs=4,
#'                               outputs=6:7)
#' result2 <- model_basic(data_example2,
#'                        orientation="io",
#'                        rts="crs")
#' efficiencies(result2)
#' targets(result2)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_basic}}

"Libraries"


#' Data: Guo and Tanaka (2001).
#'
#' Data of 5 DMUs with two symmetric triangular fuzzy inputs, Xj = (xj, alphaj), and
#' two symmetric triangular fuzzy outputs, Yj = (yj, betaj).
#' 
#' @usage data("Guo_Tanaka_2001")
#' 
#' @format Data frame with 5 rows and 9 columns. Definition of fuzzy inputs (X)
#' and fuzzy outputs (Y):
#' \describe{
#'   \item{x1}{Input 1}
#'   \item{x2}{Input 2}
#'   \item{alpha1}{spread vector Input 1}
#'   \item{alpha2}{spread vector Input 2}
#'   \item{y1}{Output 1}
#'   \item{y2}{Output 2}
#'   \item{beta1}{spread vector Output 1}
#'   \item{beta2}{spread vector Output 2}
#' }
#' 
#' @source Guo, P.; Tanaka, H. (2001). "Fuzzy DEA: A Perceptual Evaluation Method",
#' Fuzzy Sets and Systems, 119, 149–160. \doi{10.1016/S0165-0114(99)00106-2}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' data("Guo_Tanaka_2001")
#' data_example <- make_deadata_fuzzy(Guo_Tanaka_2001,
#'                                    dmus = 1,
#'                                    inputs.mL = 2:3,
#'                                    inputs.dL = 4:5,
#'                                    outputs.mL = 6:7,
#'                                    outputs.dL = 8:9)
#' result <- modelfuzzy_guotanaka(data_example,
#'                                h = seq(0, 1, by = 0.1),
#'                                orientation = "io")
#' efficiencies(result)
#'
#' @seealso \code{\link{make_deadata_fuzzy}}, \code{\link{modelfuzzy_guotanaka}},
#' \code{\link{cross_efficiency_fuzzy}}

"Guo_Tanaka_2001"


#' Data: Leon, Liern, Ruiz and Sirvent (2003).
#'
#' Data of 8 DMUs with one symmetric triangular fuzzy inputs: Xj = (xj, alphaj), and
#' one symmetric triangular fuzzy outputs: Yj = (yj, betaj).
#' 
#' @usage data("Leon2003")
#' 
#' @format Data frame with 8 rows and 5 columns. Definition of fuzzy inputs (X)
#' and fuzzy outputs (Y):
#' \describe{
#'   \item{x1}{Input 1}
#'   \item{alpha1}{spread vector Input 1}
#'   \item{y1}{Output 1}
#'   \item{beta1}{spread vector Output 1}
#' }
#' 
#' @source Leon, T.; Liern, V. Ruiz, J.; Sirvent, I. (2003). "A Possibilistic
#' Programming Approach to the Assessment of Efficiency with DEA Models", Fuzzy
#' Sets and Systems, 139, 407–419. \doi{10.1016/S0165-0114(02)00608-5}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. Replication of results in Leon et. al (2003, p. 416)
#' data("Leon2003")
#' data_example <- make_deadata_fuzzy(Leon2003,
#'                                    dmus = 1,
#'                                    inputs.mL = 2,
#'                                    inputs.dL = 3,
#'                                    outputs.mL = 4,
#'                                    outputs.dL = 5)
#' result <- modelfuzzy_possibilistic(data_example,
#'                                    h = seq(0, 1, by = 0.1),
#'                                    orientation = "io",
#'                                    rts = "vrs")
#' efficiencies(result)
#'
#' @seealso \code{\link{make_deadata_fuzzy}}, \code{\link{modelfuzzy_possibilistic}},
#' \code{\link{cross_efficiency_fuzzy}}, \code{\link{modelfuzzy_guotanaka}}

"Leon2003"


#' Data: Kao and Liu (2003).
#'
#' Data of 24 university libraries in Taiwan with one input and five outputs.
#' 
#' @usage data("Kao_Liu_2003")
#' 
#' @format Data frame with 24 rows and 11 columns. Definition of fuzzy inputs (X)
#' and fuzzy outputs (Y):
#' \describe{
#'   \item{x1 = Patronage}{It is a weighted sum of the standardized scores of faculty,
#'   graduate students, undergraduate students, and extension students in the range of 0 and 1.}
#'   \item{y1 = Collections}{Books, serials, microforms, audiovisual works, and database.}
#'   \item{y2 = Personnel}{Classified staff, unclassified staff, and student assistants.}
#'   \item{y3 = Expenditures }{Capital expenditure, operating expenditure, and
#'   special expenditure.}
#'   \item{y4 = Buildings}{Area and seats}
#'   \item{y5 = Services}{Operating hours, attendance, circulation, communication
#'   channels, range of services, amount of services, etc.}
#'   \item{beta3_l}{lower spread vector Expenditures}
#'   \item{beta3_u}{upper spread vector Expenditures}
#'   \item{beta5_l}{lower spread vector Services}
#'   \item{beta5_u}{upper spread vector Services}
#' }
#'
#' @note There are three observations that are missing: expenditures of Library 24
#' and services of Library 22 and Library 23. Kao and Liu (2000b) represent the
#' expenditures of Library 24 by the triangular fuzzy number Y = (0.11; 0.41; 1.0).
#' The services of Library 22 and Library 23 are expressed by a same triangular
#' fuzzy number Y = (0.41; 0.69; 1.0).
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @source Kao, C., Liu, S.T. (2003). “A mathematical programming approach to
#' fuzzy efficiency ranking”, International Journal of Production Economics, 85.
#' \doi{10.1016/S0925-5273(03)00026-4}
#'
#' @examples
#' # Example. Replication of results in Kao and Liu (2003, p.152)
#' data_example <- make_deadata_fuzzy(Kao_Liu_2003,
#'                                    dmus = 1,
#'                                    inputs.mL = 2,
#'                                    outputs.mL = 3:7,
#'                                    outputs.dL = c(NA, NA, 8, NA, 10),
#'                                    outputs.dR = c(NA, NA, 9, NA, 11))
#' result <- modelfuzzy_kaoliu(data_example,
#'                             kaoliu_modelname = "basic",
#'                             orientation = "oo",
#'                             rts = "vrs",
#'                             alpha = 0)
#' eff <- efficiencies(result)
#' eff
#'
#' @seealso \code{\link{make_deadata_fuzzy}}, \code{\link{model_basic}}

"Kao_Liu_2003"


#' Data: Lim and Zhu (2015).
#'
#' Data of 37 R&D project proposal relating to the Turkish iron and steel industry.
#' Authors consider one input and five outputs.
#' 
#' @usage data("Lim_Zhu_2015")
#' 
#' @format Data frame with 37 rows and 7 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1 = Budget}{Budget}
#'   \item{y1 = Indirect_economic}{Indirect economic contribution}
#'   \item{y2 = Direct_economic}{Direct economic contribution}
#'   \item{y3 = Technical}{Technical contribution}
#'   \item{y4 = Social}{Social contribution}
#'   \item{y5 = Scientific}{Scientific contribution}
#' }
#'
#' @source Lim, S.; Zhu, J. (2015). "DEA Cross-Efficiency Under Variable Returns
#' to Scale". Journal of Operational Research Society, 66(3), p. 476-487.
#' \doi{10.1057/jors.2014.13}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. Arbitrary formulation.
#' # Input-oriented model under variable returns-to-scale.
#' data("Lim_Zhu_2015")
#' data_example <- make_deadata(Lim_Zhu_2015,
#'                              dmus = 1,
#'                              ni = 1,
#'                              no = 5)
#' cross <- cross_efficiency(data_example,
#'                           epsilon = 0,
#'                           orientation = "io",
#'                           rts = "vrs",
#'                           selfapp = TRUE,
#'                           M2 = FALSE,
#'                           M3 = FALSE)
#' cross$Arbitrary$e
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_multiplier}},
#' \code{\link{cross_efficiency}}

"Lim_Zhu_2015"


#' Data: Färe, Grosskopf and Kokkelenberg (1989).
#'
#' Data of 19 coal-fired steam-electric generating plants operating in Illinois
#' in 1978. Each plant uses 3 inputs to produce 1 output.
#' 
#' @usage data("Electric_plants")
#' 
#' @format Data frame with 18 rows and 5 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1 = Labor}{Labor average annual employment}
#'   \item{x2 = Fuel}{Fuel \eqn{10^{10}} Btu}
#'   \item{x3 = Capital}{Capital MW (fixed input)}
#'   \item{y1 = Output}{Output \eqn{10^6} Kwh}
#' }
#' @source Färe, R.; Grosskopf, S.; Kokkenlenberg, E. (1989). "Measuring Plant
#' Capacity, Utilization and Technical Change: A Nonparametric Approach".
#' International Economic Review, 30(3), 655-666.
#'
#' Simar, L.; Wilson, P.W. (1998). "Sensitivity Analysis of Efficiency Scores:
#' How to Bootstrap in Nonparametric Frontier Models". Management Science, 44(1), 49-61.
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. Replication of results in Simar and Wilson (1998, p.59)
#' data("Electric_plants")
#' data_example <- make_deadata(Electric_plants,
#'                              dmus = 1,
#'                              ni = 3,
#'                              no = 1)
#' result <- model_basic(data_example,
#'                       orientation = "io",
#'                       rts = "vrs")
#' efficiencies(result)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_basic}}

"Electric_plants"


#' Data: Hua and Bian (2007).
#'
#' Data of 30 DMUs with two desirable inputs, two desirable outputs and one udesirable output.
#' 
#' @usage data("Hua_Bian_2007")
#' 
#' @format Data frame with 30 rows and 6 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1 = D-Input1}{Desirable Input 1}
#'   \item{x2 = D-Input2}{Desirable Input 2}
#'   \item{y1 = D-Output1}{Desirable Output 1}
#'   \item{y2 = D-Output2}{Desirable Output 2}
#'   \item{y3 = UD-Output1}{Undesirable Output 1}
#' }
#' 
#' @source Hua Z.; Bian Y. (2007). DEA with Undesirable Factors. In: Zhu J.,
#' Cook W.D. (eds) Modeling Data Irregularities and Structural Complexities
#' in Data Envelopment Analysis. Springer, Boston, MA.
#' \doi{10.1007/978-0-387-71607-7_6}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#'
#' @examples
#' # Example. Replication of results in Hua and Bian (2007).
#' data("Hua_Bian_2007")
#' # The third output is an undesirable output
#' data_example <- make_deadata(Hua_Bian_2007,
#'                              ni = 2,
#'                              no = 3,
#'                              ud_outputs = 3)
#'
#' # Translation parameter (vtrans_o) is set to 1500
#' result <- model_basic(data_example,
#'                       orientation = "oo",
#'                       rts = "vrs",
#'                       vtrans_o = 1500)
#' eff <- efficiencies(result)
#' 1 / eff # results M5 in Table 6-5 (p.119)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_basic}}

"Hua_Bian_2007"


#' Data: Ruggiero (2007).
#'
#' Simulated data of 35 DMUs with two inputs and one output.
#' 
#' @usage data("Ruggiero2007")
#' 
#' @format Data frame with 35 rows and 4 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1}{Input 1}
#'   \item{x2}{Input 2}
#'   \item{y1}{Output 1}
#' }
#' 
#' @source Ruggiero J. (2007). Non-Discretionary Inputs. In: Zhu J., Cook W.D.
#' (eds) Modeling Data Irregularities and Structural Complexities in Data Envelopment
#' Analysis. Springer, Boston, MA. \doi{10.1007/978-0-387-71607-7_5}
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. Replication of results in Ruggiero (2007).
#' data("Ruggiero2007")
#' # the second input is a non-discretionary input
#' datadea <- make_deadata(Ruggiero2007,
#'                         ni = 2,
#'                         no = 1,
#'                         nd_inputs = 2)
#' result <- model_basic(datadea,
#'                       orientation = "io",
#'                       rts = "crs")
#' efficiencies(result)
#' slacks(result)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_basic}}

"Ruggiero2007"


#' Data: Fried, Knox Lovell and Schmidt (1993).
#'
#' Data of 11 DMUs with two inputs and one output.
#' 
#' @usage data("Fried1993")
#' 
#' @format Data frame with 11 rows and 4 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1}{Input 1}
#'   \item{x2}{Input 2}
#'   \item{y1}{Output 1}
#' }
#' 
#' @source Ali, A.I.; Seiford, L.M. (1993). The Mathematical Programming Approach
#' to Efficiency Analysis. In Fried, H.O.; Knox Lovell, C.A.; Schmidt, S.S.(eds.),
#' The Measurement of Productive Efficiency. Techniques and Applications. New York:
#' Oxford University Press.
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. Replication of results in Ali and (1993, p.143).
#' data("Fried1993")
#' data_example <- make_deadata(Fried1993,
#'                              ni = 2,
#'                              no = 1)
#' result <- model_basic(data_example,
#'                       orientation = "oo",
#'                       rts = "vrs")
#' efficiencies(result)
#' targets(result)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_basic}}

"Fried1993"


#' Data: Coll and Blasco (2006).
#'
#' Data of six authorized dealers with two inputs and two outputs.
#' 
#' @usage data("Coll_Blasco_2006")
#' 
#' @format Data frame with 6 rows and 5 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x1 = Employees}{Number of employees}
#'   \item{x2 = Capital}{Impairment of assets}
#'   \item{y1 = Vehicles}{Number of vehicles sold}
#'   \item{y2 = Orders}{Number of orders received at the garage}
#' }
#' 
#' @source Coll-Serrano, V.; Blasco-Blasco, O. (2006). Evaluacion de la Eficiencia
#' mediante el Análisis Envolvente de Datos. Introduccion a los Modelos Básicos.
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. How to read data with deaR
#' data("Coll_Blasco_2006")
#' data_example <- make_deadata(Coll_Blasco_2006,
#'                              dmus = 1,
#'                              ni = 2,
#'                              no = 2)
#'
#' @seealso \code{\link{make_deadata}}

"Coll_Blasco_2006"


#' Data: Coelli, Rao and Battese (1998).
#'
#' Data of five DMUs with two inputs and one output. Prices for inputs are available.
#' Price for output is not from Coelli et al. (1998).
#' 
#' @usage data("Coelli_1998")
#' 
#' @format Data frame with 6 rows and 5 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{Input1}{Input 1}
#'   \item{Input2}{Input 2}
#'   \item{Output}{Output}
#'   \item{Price_input1}{Price input 1}
#'   \item{Price_input2}{Price input 2}
#'   \item{Price_output}{Price output}
#' }
#' 
#' @source Coelli, T.; Prasada Rao, D.S.; Battese, G.E. An introduction to
#' efficiency and productivity analysis. Boston: Kluwer Academic Publishers.
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. Replication of results in Coelli et al. (1998, p.166).
#' # Cost efficiency model.
#' data("Coelli_1994")
#' # Selection of prices: data_prices is the trasnpose where the prices for inputs are. 
#' data_prices <- t(Coelli_1998[, 5:6]) 
#' 
#' data_example <- make_deadata(Coelli_1998,
#'                              dmus = 1,
#'                              ni = 2,
#'                              no = 1)
#' result <- model_profit(data_example,
#'                        price_input = data_prices,
#'                        rts = "crs", 
#'                        restricted_optimal = FALSE) 
#' # notice that the option by default is restricted_optimal=TRUE
#' efficiencies(result)
#'
#' @seealso \code{\link{make_deadata}}

"Coelli_1998"


#' Data: Tone (2003).
#'
#' Data of 9 DMUs producing 2 outputs, being second output undesirable, by using 1 input.
#'  
#' @usage data("Tone2003")
#' 
#' @format Data frame with 9 rows and 4 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{x}{Input}
#'   \item{yg}{Output1 ("good" output)}
#'   \item{yb}{Output2 (undesirable "bad" output)}
#' }
#' 
#' @source Tone, K. (2003). "Dealing with undesirable outputs in DEA: A Slacks-Based
#' Measure (SBM) approach", GRIPS Research Report Series I-2003-0005.
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. Replication of results in Tone (2003), pp 10-11.
#' data("Tone2003")
#' data_example <- make_deadata(Tone2003,
#'                              ni = 1,
#'                              no = 2,
#'                              ud_outputs = 2)
#' result <- model_sbmeff(data_example,
#'                        rts = "vrs")
#' efficiencies(result)
#' targets(result)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_sbmeff}}

"Tone2003"


#' Data: Coelli, Grifell-Tatje, and Perelman (2002).
#'
#' Data of 28 airlines with 2 outputs and 4 inputs.
#'  
#' @usage data("Airlines")
#' 
#' @format Data frame with 28 rows and 7 columns. Definition of outputs (Y) and inputs (X):
#' \describe{
#'   \item{y1 = Pass}{Passenger-kilometers flown}
#'   \item{y2 = Cargo}{Freight tonne-kilometers flown}
#'   \item{x1 = Lab}{Labor (number of employees)}
#'   \item{x2 = Fuel}{Fuel (millions of gallons)}
#'   \item{x3 = Matl}{Other inputs (millions of U.S. dollar equivalent)
#'   consisting of operating and maintenance expenses excluding labor and fuel expenses}
#'   \item{x4 = Cap}{Capital (sum of the maximum takeoff weights of all aircraft
#'   flown multiplied by the number of days flown)}
#' }
#' 
#' @source Coelli, T.; Griffel-Tatje, E.; Perelman, S. (2002). "Capacity Utilization
#' and Profitability: A Decomposition of Short-Run Profit Efficiency",
#' International Journal of Production Economics 79, 261–278.
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' # Example. Replication of results in Aparicio et al. (2007).
#' data("Airlines")
#' data_example <- make_deadata(Airlines,
#'                              inputs = 4:7,
#'                              outputs = 2:3)
#' result <- model_sbmeff(data_example)
#' efficiencies(result)
#' result2 <- model_sbmeff(data_example,
#'                         kaizen = TRUE)
#' efficiencies(result2)
#'
#' @seealso \code{\link{make_deadata}}, \code{\link{model_sbmeff}}

"Airlines"


#' Data: Grifell-Tatjé and Lovell (1999).
#'
#' Data of 8 DMUs producing 1 output (Y) by using 1 input (X) for two periods of time.
#'  
#' @usage data("Grifell_Lovell_1999")
#' 
#' @format Data frame with 16 rows and 4 columns. Definition of inputs (X) and outputs (Y):
#' \describe{
#'   \item{X}{Input}
#'   \item{Y}{Output}
#' }
#' 
#' @source Grifell-Tatjé, E.; Lovel, C.A.K. (1999). "A Generalized Malmquist
#' productivity index". Top, 7(1), 81-101.  
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' 
#' # Example. Replication of results in Grifell-Tatjé and Lovell (1999, p. 100).
#' data("Grifell_Lovell_1999")
#' data_example <- make_malmquist(Grifell_Lovell_1999,
#'                                percol = 1,
#'                                dmus = 2,
#'                                inputs = 3,
#'                                outputs = 4,
#'                                arrangement = "vertical")
#'
#' result_fgnz <- malmquist_index(data_example,
#'                                orientation = "oo",
#'                                rts = "vrs",
#'                                type1 = "cont",
#'                                type2 = "fgnz")
#' 
#' result_fgnz$mi
#'                                                              
#'
#' @seealso \code{\link{make_malmquist}}, \code{\link{malmquist_index}}

"Grifell_Lovell_1999"


#' Data: Fuzzy data reading example.
#'
#' Synthetic dataset of 5 DMUs with 3 inputs and 3 outputs containing fuzzy and
#' crisp data.
#' @usage data("FuzzyExample")
#' 
#' @format Data frame with 5 rows and 15 columns.
#' \describe{
#'   \item{DMU}{DMU names}
#'   \item{Input1.mL}{First Input (crisp numbers)}
#'   \item{Input2.mL}{Second Input (left centers)}
#'   \item{Input2.mR}{Second Input (right centers)}
#'   \item{Input2.dL}{Second Input (left radii)}
#'   \item{Input2.dR}{Second Input (right radii)}
#'   \item{Input3.mL}{Third Input (centers)}
#'   \item{Input3.dL}{Third Input (radii)}
#'   \item{Output1.mL}{First Output (crisp numbers)}
#'   \item{Output2.mL}{Second Output (left centers)}
#'   \item{Output2.mR}{Second Output (right centers)}
#'   \item{Output2.dL}{Second Output (radii)}
#'   \item{Output3.mL}{Third Output (centers)}
#'   \item{Output3.dL}{Third Output (left radii)}
#'   \item{Output3.dR}{Third Output (right radii)}
#' }
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolos} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benitez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @examples
#' 
#' # Example. Reading the data.
#' data("FuzzyExample")
#' datafuzzy <- make_deadata_fuzzy(FuzzyExample, 
#'                                 inputs.mL = c(2, 3, 7),
#'                                 inputs.mR = c(NA, 4, NA),
#'                                 inputs.dL = c(NA, 5, 8),
#'                                 inputs.dR = c(NA, 6, NA),
#'                                 outputs.mL = c(9, 10 , 13),
#'                                 outputs.mR = c(NA, 11, NA),
#'                                 outputs.dL = c(NA, 12, 14),
#'                                 outputs.dR = c(NA, NA, 15))
#'                                                              
#'
#' @seealso \code{\link{make_deadata_fuzzy}}

"FuzzyExample"
