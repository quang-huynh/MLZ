#' Male Nephrops FU 28-29 (for MLeffort)
#'
#' An S4 object containing length and effort time series and life history parameters of male Nephrops in FU 28-29.
#'
#' @format An object of class \code{\linkS4class{MLZ_data}}.
#'
#' @references
#' Then, A.Y, Hoenig, J.M, and Huynh, Q.C. In revision. Estimating fishing and natural mortality rates,
#' and catchability coefficient, from a series of observations on mean length and fishing effort.
#' ICES Journal of Marine Science.
#'
#' @examples
#' data(Nephrops)
#'
"Nephrops"

#' Puerto Rico Snapper (for MLmulti)
#'
#' Mean lengths and life history for 3 species in the Puerto Rico Deepwater Snapper Complex (Unit 1):
#' silk snapper, blackfin snapper, and vermilion snapper.
#'
#' @format A list containing objects of class \code{\linkS4class{MLZ_data}}.
#'
#' @references
#' Huynh, Q.C, Gedamke, T., Hoenig, J.M, and Porch C. In press. Multispecies Extensions to a Nonequilibrium
#' Length-Based Mortality Estimator. Marine and Coastal Fisheries.
#'
#' @examples
#' data(PRSnapper)
#'
"PRSnapper"



#' Puerto Rico Mutton Snapper (for ML, MLCR)
#'
#' Mean lengths, CPUE, and life history for Puerto Rico mutton snapper.
#'
#' @format An object of class \code{\linkS4class{MLZ_data}}.
#'
#' @references
#' Huynh, Q.C., Gedamke, T., Porch, C.E., Hoenig, J.M., Walter, J.F, Bryan, M, and Brodziak, J.
#' In revision. Estimating Total Mortality Rates from Mean Lengths and Catch Rates in Non-equilibrium
#' Situations. Transactions of the American Fisheries Society.
#'
#' @examples
#' data(MuttonSnapper)
#'
"MuttonSnapper"


#' Goosefish: Northern Management Region (for ML)
#'
#' Mean lengths and life history for goosefish.
#'
#' @format An object of class \code{\linkS4class{MLZ_data}}.
#'
#' @references
#' Gedamke, T. and Hoenig, J.M. 2006. Estimating mortality from mean length data in
#' nonequilibrium situations, with application to the assessment of goosefish.
#' Transactions of the American Fisheries Society 135:476-487.
#'
#' @examples
#' data(Goosefish); Goosefish
#'
"Goosefish"


#' Silk Snapper
#'
#' Length observed from the Puerto Rico Silk Snapper handline fishery.
#'
#' @format A data frame.
#'
#' @references
#' Huynh, Q.C, Gedamke, T., Hoenig, J.M, and Porch C. In press. Multispecies Extensions to a Nonequilibrium
#' Length-Based Mortality Estimator. Marine and Coastal Fisheries.
#'
#' @examples
#' \dontrun{
#' data(SilkSnapper)
#' new("MLZ_data", Len_df = SilkSnapper)
#' }
"SilkSnapper"
