#' @title Quadratically Constrained CRS Generalized Oriented DEA model.
#'
#' @description It solves quadratically constrained CRS generalized oriented DEA
#' models, using alabama solver. By default, models are solved in a two-stage
#' process (slacks are maximized).
#'
#' @usage model_qgo(datadea,
#'             dmu_eval = NULL,
#'             dmu_ref = NULL,
#'             d_input = 1,
#'             d_output = 1,
#'             rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
#'             L = 1,
#'             U = 1,
#'             give_X = TRUE,
#'             maxslack = TRUE,
#'             weight_slack_i = 1,
#'             weight_slack_o = 1,
#'             returnqp = FALSE,
#'             ...)
#'
#' @param datadea A \code{deadata} object with \code{n} DMUs, \code{m} inputs and \code{s}
#' outputs.
#' @param dmu_eval A numeric vector containing which DMUs have to be evaluated.
#' If \code{NULL} (default), all DMUs are considered.
#' @param dmu_ref A numeric vector containing which DMUs are the evaluation
#' reference set.
#' If \code{NULL} (default), all DMUs are considered.
#' @param d_input A value, vector of length \code{m}, or matrix \code{m} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the input orientation parameters.
#' If \code{d_input} == 1 (default) and \code{d_output} == 0, it is equivalent
#' to input oriented.
#' @param d_output A value, vector of length \code{s}, or matrix \code{s} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the output orientation parameters.
#' If \code{d_input} == 0 and \code{d_output} == 1 (default), it is equivalent
#' to output oriented.
#' @param rts A string, determining the type of returns to scale, equal to "crs" (constant),
#' "vrs" (variable), "nirs" (non-increasing), "ndrs" (non-decreasing) or "grs" (generalized).
#' @param L Lower bound for the generalized returns to scale (grs).
#' @param U Upper bound for the generalized returns to scale (grs).
#' @param give_X Logical. If it is \code{TRUE}, it uses an initial vector (given by
#' the evaluated DMU) for the solver, except for "cccp". If it is \code{FALSE}, the initial vector is given
#' internally by the solver and it is usually randomly generated.
#' @param maxslack Logical. If it is \code{TRUE}, it computes the max slack solution.
#' @param weight_slack_i A value, vector of length \code{m}, or matrix \code{m} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the weights of the input slacks
#' for the max slack solution.
#' @param weight_slack_o A value, vector of length \code{s}, or matrix \code{s} x \code{ne}
#' (where \code{ne} is the length of \code{dmu_eval}) with the weights of the output
#' slacks for the max slack solution.
#' @param returnqp Logical. If it is \code{TRUE}, it returns the quadratic problems
#' (objective function and constraints) of stage 1.
#' @param ... Other parameters, like the initial vector \code{X}, to be passed to the solver.
#' 
#' @returns A list of class \code{dea} with the results for the evaluated DMUs (\code{DMU} component,
#' we note that we call "targets" to the "efficient projections"
#' in the strongly efficient frontier),
#'  along with any other necessary information to replicate the results, such as
#'  the name of the model and parameters \code{d_input}, \code{d_output}, \code{rts},
#'  \code{dmu_eval} and \code{dmu_ref}.
#'
#' @author
#' \strong{Vicente Coll-Serrano} (\email{vicente.coll@@uv.es}).
#' \emph{Quantitative Methods for Measuring Culture (MC2). Applied Economics.}
#'
#' \strong{Vicente Bolós} (\email{vicente.bolos@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Rafael Benítez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#'
#' @references
#' 
#' "A new family of models with generalized orientation in data envelopment
#' analysis". V. J. Bolós, R. Benítez, V. Coll-Serrano. International
#' Transactions in Operational Research. \doi{10.1111/itor.70063}
#'
#' @examples
#' 
#' data("PFT1981") 
#' # Selecting DMUs in Program Follow Through (PFT)
#' PFT <- PFT1981[1:49, ] 
#' PFT <- make_deadata(PFT,
#'                     inputs = 2:6, 
#'                     outputs = 7:9 )
#' eval_pft <- model_qgo(PFT, dmu_eval = 1:5)
#' efficiencies(eval_pft)
#'
#' @seealso \code{\link{model_basic}}, \code{\link{model_dir}}, \code{\link{model_lgo}}
#'
#' @import optiSolve lpSolve
#'
#' @export

model_qgo <-
  function(datadea,
           dmu_eval = NULL,
           dmu_ref = NULL,
           d_input = 1,
           d_output = 1,
           rts = c("crs", "vrs", "nirs", "ndrs", "grs"),
           L = 1,
           U = 1,
           give_X = TRUE,
           maxslack = TRUE,
           weight_slack_i = 1,
           weight_slack_o = 1,
           returnqp = FALSE,
           ...) {

    # Cheking whether datadea is of class "deadata" or not...
    if (!is.deadata(datadea)) {
      stop("Data should be of class deadata. Run make_deadata function first!")
    }

    # Checking undesirable inputs/outputs
    if (!is.null(datadea$ud_inputs) || !is.null(datadea$ud_outputs)) {
      warning("This model does not take into account the undesirable feature for inputs/outputs.")
    }

    # Checking non controllable inputs/outputs
    if (!is.null(datadea$nc_inputs) || !is.null(datadea$nc_outputs)) {
      warning("This model does not take into account the non controllable feature for inputs/outputs.")
    }

    # Checking non discretionary inputs/outputs
    if (!is.null(datadea$nd_inputs) || !is.null(datadea$nd_outputs)) {
      warning("This model does not take into account the non discretionary feature for inputs/outputs.")
    }

    # Checking rts
    rts <- tolower(rts)
    rts <- match.arg(rts)

    if (rts == "grs") {
      if (L > 1) {
        stop("L must be <= 1.")
      }
      if (U < 1) {
        stop("U must be >= 1.")
      }
    }

    # Checking solver
    solver <- "alabama"

    dmunames <- datadea$dmunames
    nd <- length(dmunames) # number of dmus

    if (is.null(dmu_eval)) {
      dmu_eval <- 1:nd
    } else if (!all(dmu_eval %in% (1:nd))) {
      stop("Invalid set of DMUs to be evaluated (dmu_eval).")
    }
    names(dmu_eval) <- dmunames[dmu_eval]
    nde <- length(dmu_eval)

    if (is.null(dmu_ref)) {
      dmu_ref <- 1:nd
    } else if (!all(dmu_ref %in% (1:nd))) {
      stop("Invalid set of reference DMUs (dmu_ref).")
    }
    names(dmu_ref) <- dmunames[dmu_ref]
    ndr <- length(dmu_ref)

    input <- datadea$input
    output <- datadea$output
    inputnames <- rownames(input)
    outputnames <- rownames(output)
    ni <- nrow(input)
    no <- nrow(output)

    namevar1 <- c(paste("lambda", 1:ndr, sep = "_"),
                  paste("phi", 1:no, sep = "_"),
                  "beta")
    namevar2 <- c(paste("lambda", 1:ndr, sep = "_"),
                  paste("slack_I", 1:ni, sep = "_"),
                  paste("slack_O", 1:no, sep = "_"))

    if (is.matrix(d_input)) {
      if ((nrow(d_input) != ni) || (ncol(d_input) != nde)) {
        stop("Invalid input orientation matrix (number of inputs x number of evaluated DMUs).")
      }
    } else if ((length(d_input) == 1) || (length(d_input) == ni)) {
      d_input <- matrix(d_input, nrow = ni, ncol = nde)
    } else {
      stop("Invalid input orientation vector.")
    }
    if (any(d_input < 0)) {
      stop("Input orientation parameters must be non-negative.")
    }
    rownames(d_input) <- inputnames
    colnames(d_input) <- dmunames[dmu_eval]

    if (is.matrix(d_output)) {
      if ((nrow(d_output) != no) || (ncol(d_output) != nde)) {
        stop("Invalid output orientation matrix (number of outputs x number of evaluated DMUs).")
      }
    } else if ((length(d_output) == 1) || (length(d_output) == no)) {
      d_output <- matrix(d_output, nrow = no, ncol = nde)
    } else {
      stop("Invalid output orientation vector.")
    }
    if (any(d_output < 0)) {
      stop("Output orientation parameters must be non-negative.")
    }
    rownames(d_output) <- outputnames
    colnames(d_output) <- dmunames[dmu_eval]

    inputref <- matrix(input[, dmu_ref], nrow = ni)
    outputref <- matrix(output[, dmu_ref], nrow = no)

    target_input <- NULL
    target_output <- NULL
    orientation_param <- NULL

    DMU <- vector(mode = "list", length = nde)
    names(DMU) <- dmunames[dmu_eval]

    ###########################

    # Objective function coefficients stage 1
    f.obj <- linfun(a = c(rep(0, ndr + no), 1), id = namevar1)

    # Lower and upper bounds constraints stage 1
    lbcon1 <- lbcon(val = 0, id = namevar1)
    ubcon1 <- NULL

    if (rts == "crs") {
      f.con.rs <- NULL # Stage 1
      f.con2.rs <- NULL # Stage 2
      f.dir.rs <- NULL
      f.rhs.rs <- NULL
    } else {
      f.con.rs <- cbind(matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = no + 1))
      f.con2.rs <- cbind(matrix(1, nrow = 1, ncol = ndr), matrix(0, nrow = 1, ncol = ni + no))
      f.rhs.rs <- 1
      if (rts == "vrs") {
        f.dir.rs <- "="
        ubcon1 <- ubcon(val = rep(1, ndr), id = namevar1[1:ndr])
      } else if (rts == "nirs") {
        f.dir.rs <- "<="
        ubcon1 <- ubcon(val = rep(1, ndr), id = namevar1[1:ndr])
      } else if (rts == "ndrs") {
        f.dir.rs <- ">="
      } else {
        f.con.rs <- rbind(f.con.rs, f.con.rs)
        f.con2.rs <- rbind(f.con2.rs, f.con2.rs)
        f.dir.rs <- c(">=", "<=")
        f.rhs.rs <- c(L, U)
        ubcon1 <- ubcon(val = rep(U, ndr), id = namevar1[1:ndr])
      }
    }

    if (maxslack && (!returnqp)) {

      # Checking weights
      if (is.matrix(weight_slack_i)) {
        if ((nrow(weight_slack_i) != ni) || (ncol(weight_slack_i) != nde)) {
          stop("Invalid weight input matrix (number of inputs x number of evaluated DMUs).")
        }
      } else if ((length(weight_slack_i) == 1) || (length(weight_slack_i) == ni)) {
        weight_slack_i <- matrix(weight_slack_i, nrow = ni, ncol = nde)
      } else {
        stop("Invalid weight input vector (number of inputs).")
      }
      rownames(weight_slack_i) <- inputnames
      colnames(weight_slack_i) <- dmunames[dmu_eval]

      if (is.matrix(weight_slack_o)) {
        if ((nrow(weight_slack_o) != no) || (ncol(weight_slack_o) != nde)) {
          stop("Invalid weight output matrix (number of outputs x number of evaluated DMUs).")
        }
      } else if ((length(weight_slack_o) == 1) || (length(weight_slack_o) == no)) {
        weight_slack_o <- matrix(weight_slack_o, nrow = no, ncol = nde)
      } else {
        stop("Invalid weight output vector (number of outputs).")
      }
      rownames(weight_slack_o) <- outputnames
      colnames(weight_slack_o) <- dmunames[dmu_eval]

      # Constraints matrix stage 2
      f.con2.1 <- cbind(inputref, diag(ni), matrix(0, nrow = ni, ncol = no))
      f.con2.2 <- cbind(outputref, matrix(0, nrow = no, ncol = ni), -diag(no))
      f.con2 <- rbind(f.con2.1, f.con2.2, f.con2.rs)

      # Directions vector stage 2
      f.dir2 <- c(rep("=", ni + no), f.dir.rs)

    }

    # Linear Directions vector stage 1
    f.dir <- c(rep("<=", ni), rep(">=", no), f.dir.rs)

    for (i in 1:nde) {

      ii <- dmu_eval[i]

      # Linear Constraints matrix stage 1
      f.con.1 <- cbind(inputref, matrix(0, nrow = ni, ncol = no), input[, ii] * d_input[, i])
      f.con.2 <- cbind(outputref, -output[, ii] * diag(no), matrix(0, nrow = no, ncol = 1))
      f.con <- rbind(f.con.1, f.con.2, f.con.rs)

      # Linear Right hand side vector stage 1
      f.rhs <- c(input[, ii], rep(0, no), f.rhs.rs)
      
      rownames(f.con) <- paste("lc", 1:nrow(f.con), sep = "") # to prevent names errors in lincon
      lincon1 <- lincon(A = f.con, dir = f.dir, val = f.rhs, id = namevar1)
      
      # Quadratic constraints stage 1
      qclist <- vector("list", no)
      for (ir in 1:no) {
        qcmat <- matrix(0, nrow = ndr + no + 1, ncol = ndr + no + 1)
        qcmat[ndr + no + 1, ndr + ir] <- d_output[ir, i] / 2
        qcmat[ndr + ir, ndr + no + 1] <- d_output[ir, i] / 2
        avec <- rep(0, ndr + no + 1)
        avec[ndr + ir] <- -1
        qclist[[ir]] <- quadcon(Q = qcmat, a = avec, val = -1, id = namevar1)
      }
      names(qclist) <- paste("qc", 1:no, sep = "")

      mycop <- cop(f = f.obj, max = TRUE, lb = lbcon1, ub = ubcon1, lc = lincon1)
      mycop$qc <- qclist

      if (returnqp) {

        DMU[[i]] <- mycop

      } else {

        # Initial vector
        if ((ii %in% dmu_ref) && give_X) {
          Xini <- c(rep(0, ndr), rep(1, no), 1)
          Xini[which(dmu_ref == ii)] <- 1
          names(Xini) <- namevar1
        } else {
          Xini <- NULL
        }

        res <- solvecop(op = mycop, solver = solver, quiet = TRUE, X = Xini, ...)

        if (res$status == "successful convergence") {

          res <- res$x
          beta <- res[ndr + no + 1]

          rho <- (1 - sum(d_input[, i]) * beta / ni) * no /
            (sum(1/(1 - beta * d_output[, i])))
          names(rho) <- "rho"

          proj_input <- input[, ii] * (1 - beta * d_input[, i])
          names(proj_input) <- inputnames
          proj_output <- output[, ii] / (1 - beta * d_output[, i])
          names(proj_output) <- outputnames

          if (maxslack) {

            # Objective function coefficients stage 2
            f.obj2 <- c(rep(0, ndr), weight_slack_i[, i], weight_slack_o[, i])

            # Right hand side vector stage 2
            f.rhs2 <- c(proj_input, proj_output, f.rhs.rs)

            res <- lp("max", f.obj2, f.con2, f.dir2, f.rhs2)$solution

          }

          lambda <- res[1 : ndr]
          names(lambda) <- dmunames[dmu_ref]

          target_input <- as.vector(inputref %*% lambda)
          names(target_input) <- inputnames
          target_output <- as.vector(outputref %*% lambda)
          names(target_output) <- outputnames

          slack_input <- proj_input - target_input
          names(slack_input) <- inputnames
          slack_output <- target_output - proj_output
          names(slack_output) <- outputnames
          
          # Cambio de notación: El target es lo que era la proyección y la proyección
          # eficiente es lo que era el target.
          effproj_input <- target_input
          effproj_output <- target_output
          target_input <- proj_input
          target_output <- proj_output

        } else {
          
          rho <- NA
          beta <- NA
          lambda <- NA
          target_input <- NA
          target_output <- NA
          slack_input <- NA
          slack_output <- NA
          effproj_input <- NA
          effproj_output <- NA
          
        }

        DMU[[i]] <- list(efficiency = rho,
                         beta = beta,
                         lambda = lambda,
                         target_input = target_input, target_output = target_output,
                         slack_input = slack_input, slack_output = slack_output,
                         effproj_input = effproj_input, effproj_output = effproj_output)

      }

    }

    deaOutput <- list(modelname = "qgo",
                      d_input = d_input,
                      d_output = d_output,
                      rts = rts,
                      L = L,
                      U = U,
                      DMU = DMU,
                      data = datadea,
                      rho = rho,
                      beta = beta,
                      dmu_eval = dmu_eval,
                      dmu_ref = dmu_ref,
                      maxslack = maxslack,
                      weight_slack_i = weight_slack_i,
                      weight_slack_o = weight_slack_o)

    return(structure(deaOutput, class = "dea"))

  }
