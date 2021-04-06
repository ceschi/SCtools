#### Functions

gen_placebos_l1o <- function(dataprep.out,
															synth.out,
															l1o = T,
															#l1o = T just leave treated out of the placebo
															#l1o = F og treated is blended in the placebo
															Sigf.ipop = 5,
															strategy = "sequential",
															.optimethod = "Nelder-Mead") {
	
	# Inputs Checkss
	if(!all(names(dataprep.out) %in% dataprep_object_names)){
		stop("Invalid dataprep object. Have you run `dataprep` on your data?")
	}
	
	if(!all(names(synth.out) %in% synth_object_names)){
		stop("Invalid synth output. Have you run the `synth` function?")
	}
	
	if(is.integer(Sigf.ipop) | Sigf.ipop <= 0){
		stop("You have not passed a valid argument for Signf.ipop. Please pass a positive integer")
	}
	
	strategy_match <- match.arg(strategy, c("sequential", "multiprocess", "multicore", "multisession"))
	
	unit.numbers <- NULL
	
	tr <- as.numeric(dataprep.out$tag$treatment.identifier)
	
	# plain 
	if (l1o){
		n <- length(dataprep.out$tag$controls.identifier)
		names.and.numbers <-
				subset(dataprep.out$names.and.numbers, unit.numbers != tr)
	} 
	
	# this one keeps also the original treated
	if (!l1o){
		n <- length(dataprep.out$tag$controls.identifier) + 1
		names.and.numbers <- dataprep.out$names.and.numbers
	} 
	
	b <-
		data.frame(matrix(
			0,
			ncol = n,
			nrow = length(dataprep.out$tag$time.plot)
		))
	
	strategy <- paste0("future::", strategy_match)
	
	plan(strategy)
	oplan <- plan()
	
	mspe2 <- furrr::future_map(1:n, ~syn_plac_l1o(.x, 
																						dataprep.out, 
																						Sigf.ipop, 
																						tr, names.and.numbers,
																						optimethod = .optimethod,
																						l1o = l1o), 
														 .progress = T)
	
	b <- dplyr::bind_cols(purrr::map(mspe2,"a"))
	b <- setNames(b, paste('synthetic', as.character(names.and.numbers[ ,2]), sep = '.'))
	
	mspe.placs <- purrr::map(mspe2, "s.mspe")
	mspe.placs <- as.data.frame(unlist(mspe.placs))
	
	# ww <- dplyr::bind_rows(purrr::map(mspe2, "weights"))
	# names(ww) <- paste0('plac.weights.', names.and.numbers[ ,2])
	# row.names(ww) <- names.and.numbers[,2]
	# ww <- reshape(idvar = 'pool',
	#               timevar = 'treat',
	#               direction = 'wide')
	# ww <- ww[order(ww$pool),]
	# ww$pool <- NULL
	ww <- dplyr::bind_cols(purrr::map(mspe2, 'weights'))
	
	on.exit(plan(oplan), add = TRUE)
	
	df <-
		cbind(
			b,
			dataprep.out$Y0,
			dataprep.out$Y1,
			dataprep.out$Y0plot %*% synth.out$solution.w,
			dataprep.out$tag$time.plot
		)
	colnames(df)[(ncol(df) - 2):ncol(df)] <-c('Y1', 'synthetic.Y1', 'year')
	
	t0 <- as.numeric(dataprep.out$tag$time.optimize.ssr[1])
	t1 <-
		as.numeric(dataprep.out$tag$time.optimize.ssr[length(dataprep.out$tag$time.optimize.ssr)]) +
		1
	treated.name <-
		as.character(dataprep.out$names.and.numbers$unit.names[dataprep.out$names.and.numbers[, 2] %in% dataprep.out$tag$treatment.identifier])
	
	loss.v <- synth.out$loss.v
	
	res2 <-
		list(
			df = df,
			weights.w = as.matrix(ww),
			mspe.placs = mspe.placs,
			t0 = t0,
			t1 = t1,
			tr = tr,
			names.and.numbers = names.and.numbers,
			n = n,
			treated.name = treated.name,
			loss.v = loss.v
		)
	
	class(res2) <- append(class(res2),"tdf")
	
	return(res2)
}


syn_plac_l1o <- function(i, 
										 dataprep.out, 
										 Sigf.ipop, 
										 tr, 
										 names.and.numbers, 
										 optimethod = "Nelder-Mead",
										 l1o = T) {
	
	nunits <- length(dataprep.out$tag$foo)
	
	if (!l1o){
		# bind treated (1) with controls (0)
		X <- cbind(dataprep.out$X1, dataprep.out$X0)
		Z <- cbind(dataprep.out$Z1, dataprep.out$Z0)
		
		X0 <- X[, -i]
		X1 <- matrix(X[, i, drop = F])
		
		Z0 <- Z[, -i]
		Z1 <- matrix(Z[, i, drop = F])
		
		Yplot <- cbind(dataprep.out$Y1plot, dataprep.out$Y0plot)
		
		Y0plot <- Yplot[, -i]
		Y1plot <- matrix(Yplot[, i, drop = F])
		
		nunits <- nunits + 1 
		
		foo <- dataprep.out$tag$foo
		
		controls.id <- c(dataprep.out$tag$treatment.identifier,
										 dataprep.out$tag$controls.identifier)[-i]
		treated.id <- c(dataprep.out$tag$treatment.identifier,
										dataprep.out$tag$controls.identifier)[i]
		
	} else {
		
		X0 <- dataprep.out$X0[, -i]
		X1 <- matrix(dataprep.out$X0[, i, drop = F])
		Z0 <- dataprep.out$Z0[, -i]
		Z1 <- matrix(dataprep.out$Z0[, i, drop = F])
		Y0plot <- dataprep.out$Y0plot[, -i]
		Y1plot <- dataprep.out$Y0plot[, i, drop = F]
		
		# foo is the original panel, oddly stored as char
		# create empty one to remove from it the truly treated unit
		foo <- character(length = nunits)
		
		# where unit identifiers are
		id <- as.numeric(dataprep.out$tag$unit.variable)
		
		# finds in foo the treated block coming from the main SCM
		d <-
			which(regexpr(tr, stringr::str_split(dataprep.out$tag$foo[id],
																					 ', ')[[1]]) == T)
		
		# fill up the empty foo without obs from treated unit
		for (j in 1:length(dataprep.out$tag$foo)) {
			foo[j] <-
				paste(stringr::str_split(dataprep.out$tag$foo[j],
																 ', ')[[1]][-c(d[1]:d[length(d)])], 
							collapse =', ')
		}
			controls.id <- dataprep.out$tag$controls.identifier[-i]
			treated.id <- dataprep.out$tag$controls.identifier[i]
	}
	
	#TODO: reassemble properly datapreps and foo objs
	
	# recycle old info
	predictors <- dataprep.out$tag$predictors
	predictors.op <- dataprep.out$tag$predictors.op
	special.predictors <- dataprep.out$tag$special.predictors
	dependent <- dataprep.out$tag$dependent
	unit.variable <- dataprep.out$tag$unit.variable
	time.variable <- dataprep.out$tag$time.variable
	time.predictors.prior <- dataprep.out$tag$time.predictors.prior
	time.optimize.ssr <- dataprep.out$tag$time.optimize.ssr
	time.plot <- dataprep.out$tag$time.plot
	unit.names.variable <- dataprep.out$tag$unit.names.variable
	
	tag <- list(
		foo = as.character(foo),
		predictors = predictors,
		predictors.op = predictors.op,
		special.predictors = special.predictors,
		dependent = dependent,
		unit.variable = unit.variable,
		time.variable = time.variable,
		# next two are mod'd
		treatment.identifier = treated.id,
		controls.identifier = controls.id,
		time.predictors.prior = time.predictors.prior,
		time.optimize.ssr = time.optimize.ssr,
		time.plot = time.plot,
		unit.names.variable = unit.names.variable
	)
	dp <- list(
		X0 = X0,
		X1 = X1,
		Z0 = Z0,
		Z1 = Z1,
		Y0plot = Y0plot,
		Y1plot = Y1plot,
		names.and.numbers = names.and.numbers,
		tag = tag
	)
	s.out <- Synth::synth(data.prep.obj = dp, 
												Sigf.ipop = Sigf.ipop, optimxmethod = optimethod)
	
	# a is the synth for i-th unit in the placebo
	a <- data.frame(dp$Y0plot %*% s.out$solution.w)
	
	# no big fuss on it
	# wei <- data.frame(weights = s.out$solution.w,
	# 									pool = as.numeric(row.names(s.out$solution.w)),
	# 									treat = cases[i],
	# 									row.names = NULL)
	
	wei <- s.out$solution.w
	
	s.mspe <- s.out$loss.v
	res <- list(a = a, 
							s.mspe = s.mspe, 
							weights = wei)
	return(res)
}
