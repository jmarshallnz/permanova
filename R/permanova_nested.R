
# functions for nested ANOVA

get_ss_nested <- function(diss, A, B_within_A)
{
  n_tot    <- length(B_within_A)
  
  # compute SS_within_A
  A_vals <- unique(A)
  
  SS_within_A <- 0
  for (a in A_vals)
  {
    wch <- A == a
    SS_within_A <- SS_within_A + sum(diss[wch, wch]) / (2*sum(wch))
  }
  
  # now compute SS_res
  b_in_a <- paste(A, B_within_A)
  SS_res <- 0
  B_vals <- unique(b_in_a)
  for (b in B_vals)
  {
    wch <- b_in_a == b
    SS_res <- SS_res + sum(diss[wch, wch]) / (2*sum(wch))
  }
  SS_B_in_A <- SS_within_A - SS_res
  
  return(list(SS_within_A = SS_within_A,
              SS_B_in_A = SS_B_in_A,
              SS_res = SS_res))
}

compute_F_nested <- function(diss, A, B_within_A)
{
  n_tot    <- length(B_within_A)
  
  # compute SS_within_A
  A_vals <- unique(A)
  
  SS_within_A <- 0
  for (a in A_vals)
  {
    wch <- A == a
    SS_within_A <- SS_within_A + sum(diss[wch, wch]) / (2*sum(wch))
  }
  
  # now compute SS_res
  b_in_a <- paste(A, B_within_A)
  SS_res <- 0
  B_vals <- unique(b_in_a)
  for (b in B_vals)
  {
    wch <- b_in_a == b
    SS_res <- SS_res + sum(diss[wch, wch]) / (2*sum(wch))
  }
  SS_B_in_A <- SS_within_A - SS_res
  
  # compute F-value
  F <- SS_B_in_A / (length(B_vals) - length(A_vals)) / (SS_res / (n_tot-length(B_vals)))
  return(F)
}

get_ss_main <- function(diss, A, B_within_A) {
  n_tot  = length(B_within_A)
  A_vals = unique(A)
  b_in_a = paste(A, B_within_A)
  B_vals = unique(b_in_a)
  
  # compute SS_A
  SS_tot = sum(diss) / (2*n_tot)
  
  # compute SS_within_A
  SS_within_A <- 0
  for (a in A_vals)
  {
    wch <- A == a
    SS_within_A <- SS_within_A + sum(diss[wch, wch]) / (2*sum(wch))
  }
  # compute SS_res
  SS_res <- 0
  for (b in B_vals)
  {
    wch <- b_in_a == b
    SS_res <- SS_res + sum(diss[wch, wch]) / (2*sum(wch))
  }
  
  SS_A      = SS_tot - SS_within_A
  SS_B_in_A = SS_within_A - SS_res
  
  return(list(SS_tot = SS_tot,
              SS_A = SS_tot - SS_within_A,
              SS_within_A = SS_within_A,
              SS_B_in_A = SS_B_in_A,
              SS_res = SS_res))
}

compute_F_main <- function(diss, A, B_within_A, randomA=FALSE, randomB=TRUE, SStype=1)
{
  # TODO: Figure out what happens when
  #        A is random
  
  n_tot  = length(B_within_A)
  A_vals = unique(A)
  b_in_a = paste(A, B_within_A)
  B_vals = unique(b_in_a)
  
  # compute SS_A
  SS_tot = sum(diss) / (2*n_tot)
  
  # compute SS_within_A
  SS_within_A <- 0
  for (a in A_vals)
  {
    wch <- A == a
    SS_within_A <- SS_within_A + sum(diss[wch, wch]) / (2*sum(wch))
  }
  # compute SS_res
  SS_res <- 0
  for (b in B_vals)
  {
    wch <- b_in_a == b
    SS_res <- SS_res + sum(diss[wch, wch]) / (2*sum(wch))
  }
  
  SS_A      = SS_tot - SS_within_A
  SS_B_in_A = SS_within_A - SS_res
  SS_tot - SS_res
  # compute F-values
  MS_A      = SS_A / (length(A_vals)-1)
  MS_B_in_A = SS_B_in_A / (length(B_vals) - length(A_vals))
  MS_res    = SS_res / (n_tot - length(B_vals))
  
  if (randomB) {
    # See Sanni + Ukaegbu 2012
    
    # EMS(A)      = V(res) + k_a S(Ri) + k_b V(B_in_A)
    # EMS(B_in_A) = V(res) + k_1 V(B_in_A)
    
    # where
    n_i = table(A)
    k_a = 1/(length(n_i)-1) * (sum(n_i) - 1/sum(n_i)*sum(n_i^2))
    
    n_ij = table(A, B_within_A)
    k_b = (sum(1/n_i * rowSums(n_ij^2)) - 1/sum(n_i)*sum(n_ij^2))/(length(n_i) - 1)
    
    k_1 = (sum(n_i) - sum( 1/n_i * rowSums(n_ij^2))) / (length(B_vals)-length(A_vals))
    
    # when S(R_i) = 0 we have EMS(A) = V(res) + k_b V(B_in_A)
    # so a linear combination of these that are will involve adding on some V(res)
    # so that they're equal.
    
    num   = (k_b/k_1 - 1) * MS_res + MS_A
    denom = k_b/k_1 * MS_B_in_A
    
    F = num/denom
  } else {
    F = MS_A / MS_res
  }
  return(F)
}

# functions for permutation
get_perms_nested <- function(A, B_within_A, nperm=10000) {
  # unique values of A
  A_vals <- unique(A)
  
  b_perms_in_a <- matrix('', nrow=length(B_within_A), ncol=nperm)
  for (i in 1:nperm)
  {
    # firstly, permute items within the main effect, between the nested random effect
    b <- numeric(length(B_within_A))
    for (a in A_vals)
    {
      wch <- A == a
      b[wch] <- sample(B_within_A[wch])
    }
    b_perms_in_a[,i] <- b
  }
  b_perms_in_a
}

calc_f_nested <- function(diss, A, B_within_A, b_perms_in_a) {
  # value for our data
  F <- compute_F_nested(diss, A, B_within_A)
  
  # permutations
  nperm <- ncol(b_perms_in_a)
  Fp <- numeric(nperm)
  
  for (i in 1:ncol(b_perms_in_a)) {
    # compute the new F's
    Fp[i] <- compute_F_nested(diss, A, b_perms_in_a[,i])
  }
  list(F=F, P=sum(Fp > F)/nperm)
}

perm_nested <- function(diss, A, B_within_A, nperm=10000) {
  perms <- get_perms_nested(A, B_within_A, nperm)
  calc_f_nested(diss, A, B_within_A, perms)
}

# Main effect
get_perms_main <- function(A, B_within_A, nperm=10000) {
  
  b_within_a = paste0(A,'_',B_within_A)
  b_in_a_map = unique(cbind(b_within_a, A))
  b_val = unique(b_within_a)
  
  perms <- matrix("", nrow=length(A), ncol=nperm)
  for (i in 1:nperm)
  {
    # permute the A labels among the B's
    a_samp = sample(b_in_a_map[,2])
    
    # now assign these to each b_within_a
    a <- numeric(length(A))
    for (j in seq_len(nrow(b_in_a_map))) {
      wch = b_within_a == b_in_a_map[j,1]
      a[wch] = a_samp[j]
    }
    perms[,i] <- a
  }
  perms
}

calc_f_main <- function(diss, A, B_within_A, a_perms) {
  # value for our data
  F <- compute_F_main(diss, A, B_within_A)
  
  # permutations
  nperm <- ncol(a_perms)
  Fp <- numeric(nperm)
  
  for (i in 1:ncol(a_perms)) {
    # compute the new F's
    Fp[i] <- compute_F_main(diss, a_perms[,i], B_within_A)
  }
  list(F=F, P=sum(Fp >= F)/nperm)
}

perm_main = function(diss, A, B_within_A, nperm=10000) {
  perms <- get_perms_main(A, B_within_A, nperm)
  calc_f_main(diss, A, B_within_A, perms)
}
