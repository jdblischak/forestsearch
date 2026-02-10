
#' Cross-Validation Subgroup Match Summary
#'
#' Summarizes the match between cross-validation subgroups and analysis subgroups.
#'
#' @param sg1 Character vector. Subgroup 1 labels for each fold.
#' @param sg2 Character vector. Subgroup 2 labels for each fold.
#' @param confs Character vector. Confounder names.
#' @param sg_analysis Character vector. Subgroup analysis labels.
#'
#' @return List with indicators for any match, exact match, one match, and covariate-specific matches.
#'
#' @importFrom stringr str_sub str_length
#' @export

CV_sgs <- function(sg1, sg2, confs, sg_analysis){
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' required.")
  }

  any_found <- ifelse(!is.na(sg1) | !is.na(sg2),1,0)
  sg_depth <- length(sg_analysis)

  if(sg_depth==2){
    sg1a <- sg_analysis[1]
    sg2a <- sg_analysis[2]
    ## Exact match on both to analysis data
    exact_match <- ifelse((sg1 == sg1a | sg2 == sg1a) & (sg1 == sg2a | sg2 == sg2a),1,0)
    exact_match[is.na(exact_match)] <- 0.0
    ## At least 1 exact match
    one_match <- ifelse((sg1 == sg1a | sg2 == sg1a) | (sg1 == sg2a | sg2 == sg2a),1,0)
    one_match[is.na(one_match)] <- 0.0
    ## Cov 1 exact
    cov1_match <- ifelse((sg1 == sg1a | sg2 == sg1a),1,0)
    cov1_match[is.na(cov1_match)] <- 0.0
    cov2_match <- ifelse((sg1 == sg2a | sg2 == sg2a),1,0)
    cov2_match[is.na(cov2_match)] <- 0.0
    ## Find confounder names involved in sg1a and sg2a
    ## Add { or !{ to names for matching (a bit tedious, but let's see)
    dda <- charmatch("{",sg1a, nomatch=0)
    ddb <- charmatch("!{",sg1a, nomatch=0)
    if(dda ==1) aa <- rep("{",length(confs))
    if(ddb ==1) aa <- rep("!{",length(confs))

    temp <- paste0(aa,confs)
    temp <-paste0(temp,"}")
    ## Is sg1a confounder (NOT necessarily same cut) involved in any
    loc_name <- charmatch(temp,sg1a)
    index_name <- which(loc_name==1)

    if(length(index_name)==0){
      temp <- paste0(aa,confs)
      loc_name <- charmatch(temp,sg1a)
      index_name <- which(loc_name==1)
    }

    # Find exact match
    if(length(index_name)>1){
      confs2 <- confs[index_name]
      lc <- stringr::str_length(confs2)
      if(dda==1) ctoget <- stringr::str_sub(sg1a,2,max(lc))
      if(ddb==1) ctoget <- stringr::str_sub(sg1a,3,max(lc))
      itoget <- which(confs == ctoget)
      cfs <- confs[itoget]
    }
    if(length(index_name)==1) cfs <- confs[which(loc_name==1)]

    bb1 <- grepl(cfs,sg1)
    bb2 <- grepl(cfs,sg2)
    cov1_any <- ifelse(bb1 | bb2, 1,0)
    rm("bb1","bb2","cfs")


    ## Second
    dda <- charmatch("{",sg2a, nomatch=0)
    ddb <- charmatch("!{",sg2a, nomatch=0)
    if(dda ==1) aa <- rep("{",length(confs))
    if(ddb ==1) aa <- rep("!{",length(confs))
    temp <- paste0(aa,confs)
    temp <-paste0(temp,"}")

    loc_name <- charmatch(temp,sg2a)
    index_name <- which(loc_name==1)

    if(length(index_name)==0){
      temp <- paste0(aa,confs)
      loc_name <- charmatch(temp,sg2a)
      index_name <- which(loc_name==1)
    }

    # Find exact match
    if(length(index_name)>1){
      confs2 <- confs[index_name]
      lc <- stringr::str_length(confs2)
      if(dda==1) ctoget <- stringr::str_sub(sg2a,2,max(lc))
      if(ddb==1) ctoget <- stringr::str_sub(sg2a,3,max(lc))
      itoget <- which(confs == ctoget)
      cfs <- confs[itoget]
    }
    if(length(index_name)==1) cfs <- confs[which(loc_name==1)]

    bb1 <- grepl(cfs,sg1)
    bb2 <- grepl(cfs,sg2)
    cov2_any <- ifelse(bb1 | bb2, 1,0)
  }
  if(sg_depth==1){
    sg1a <- sg_analysis[1]
    ## Exact match on both to analysis data
    exact_match <- ifelse((sg1 == sg1a | sg2 == sg1a),1,0)
    exact_match[is.na(exact_match)] <- 0.0
    dda <- charmatch("{",sg1a, nomatch=0)
    ddb <- charmatch("!{",sg1a, nomatch=0)
    if(dda ==1) aa <- rep("{",length(confs))
    if(ddb ==1) aa <- rep("!{",length(confs))
    temp <- paste0(aa,confs)
    temp <-paste0(temp,"}")
    ## Is sg1a confounder (NOT necessarily same cut) involved in any
    loc_name <- charmatch(temp,sg1a)
    index_name <- which(loc_name==1)

    if(length(index_name)==0){
      temp <- paste0(aa,confs)
      loc_name <- charmatch(temp,sg1a)
      index_name <- which(loc_name==1)
    }

    # If names have numbers which may not be unique
    # E.g., "z1", "z11", this will match both
    # Find exact
    # Find exact match
    if(length(index_name)>1){
      confs2 <- confs[index_name]
      lc <- stringr::str_length(confs2)
      if(dda==1) ctoget <- stringr::str_sub(sg1a,2,max(lc))
      if(ddb==1) ctoget <- stringr::str_sub(sg1a,3,max(lc))
      itoget <- which(confs == ctoget)
      cfs <- confs[itoget]
    }
    if(length(index_name)==1) cfs <- confs[which(loc_name==1)]
    bb1 <- grepl(cfs,sg1)
    bb2 <- grepl(cfs,sg2)
    cov1_any <- ifelse(bb1 | bb2, 1,0)
    one_match <- exact_match
    cov2_any <- NA

    cov1_match <- exact_match
    cov2_match <- NA

  }
  return(list(any_found=any_found, exact_match=exact_match, one_match=one_match, cov1_any=cov1_any, cov2_any=cov2_any, cov1_exact=cov1_match, cov2_exact=cov2_match))
}

