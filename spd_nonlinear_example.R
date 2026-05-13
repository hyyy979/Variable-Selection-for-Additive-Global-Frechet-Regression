library(MASS)
library(Matrix)

# ============================================================
# Cholesky-based SPD utility functions
# ============================================================

# Cholesky-based inner product
cholesky_inner_product <- function(Y_i, Y_j) {
  
  L_i <- chol(Y_i)
  L_j <- chol(Y_j)
  
  return(sum((L_i - L_j)^2))
}

# Cholesky-based distance to zero
cholesky_to_zero <- function(Y_i) {
  
  L_i <- chol(Y_i)
  
  return(sum(L_i^2))
}


# ============================================================
# SCAD derivative
# ============================================================

# Right derivative of the SCAD penalty
scadrightderv <- function(tt, a = 3.7, lambda) {
  
  lambda * (tt <= lambda) +
    pmax(a * lambda - tt, 0) * (tt > lambda) / (a - 1)
}


# ============================================================
# Covariate generation
# ============================================================

# Generate correlated covariates with AR(1) structure
generate_X <- function(n, p, rho) {
  
  Sigma <- outer(
    1:p,
    1:p,
    function(i, j) rho^abs(i - j)
  )
  
  Z <- mvrnorm(
    n,
    mu = rep(0, p),
    Sigma = Sigma
  )
  
  X <- 2 * pnorm(Z)
  
  return(X)
}


# ============================================================
# Nonlinear SPD response generation
# ============================================================

# Generate nonlinear SPD-valued responses
generate_Y_spd_nonlinear <- function(
    X,
    active_mean = c(1, 3),
    active_scale = c(5, 7, 9),
    mu_0 = 3,
    beta = 2,
    sigma_0 = 1,
    gamma = 4,
    v_1 = 1,
    v_2 = 0.5
) {
  
  n <- nrow(X)
  
  X1 <- X[, active_mean[1]]
  X3 <- X[, active_mean[2]]
  
  X5 <- X[, active_scale[1]]
  X7 <- X[, active_scale[2]]
  X9 <- X[, active_scale[3]]
  
  # Nonlinear mean structure
  mean_signal <- 3 * X1^2 +
    sin(2 * pi * X3)
  
  mu_X <- rnorm(
    n,
    mean = mu_0 + beta * mean_signal,
    sd = v_1
  )
  
  # Nonlinear scale structure
  scale_signal <- exp(-X5) +
    2 * exp(-2 * (X7 - 1)^2) +
    2 / (1 + abs(X9))
  
  shape_X <- (
    sigma_0 + gamma * scale_signal
  )^2 / v_2
  
  scale_X <- v_2 / (
    sigma_0 + gamma * scale_signal
  )
  
  sigma_X <- rgamma(
    n,
    shape = shape_X,
    scale = scale_X
  )
  
  # SPD matrix construction
  I_matrix <- diag(3)
  
  U_matrix <- matrix(0, 3, 3)
  U_matrix[upper.tri(U_matrix)] <- 1
  
  Y_list <- lapply(
    1:n,
    function(i) {
      
      A <- (mu_X[i] + sigma_X[i]) * I_matrix +
        sigma_X[i] * U_matrix
      
      return(t(A) %*% A)
    }
  )
  
  return(Y_list)
}


# ============================================================
# Kernel centering
# ============================================================

# Center Gram matrix
center_gram <- function(K) {
  
  n <- nrow(K)
  
  H <- diag(n) - matrix(1, n, n) / n
  
  return(H %*% K %*% H)
}


# ============================================================
# Gaussian kernel Gram matrices
# ============================================================

# Compute Gaussian kernel Gram matrices
compute_gram_matrices <- function(X, gamma = 1) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  G_list <- vector("list", p)
  
  for (j in 1:p) {
    
    xj <- X[, j, drop = FALSE]
    
    D <- as.matrix(dist(xj))^2
    
    K <- exp(-gamma * D)
    
    G_list[[j]] <- center_gram(K)
  }
  
  return(G_list)
}


# ============================================================
# Validation Gram matrices
# ============================================================

# Compute validation kernel matrices
compute_validation_gram <- function(
    X_valid,
    X_train,
    gamma = 1
) {
  
  n_train <- nrow(X_train)
  n_valid <- nrow(X_valid)
  
  p <- ncol(X_train)
  
  G_valid_list <- vector("list", p)
  
  H_train <- diag(n_train) -
    matrix(1, n_train, n_train) / n_train
  
  H_valid <- diag(n_valid) -
    matrix(1, n_valid, n_valid) / n_valid
  
  for (j in 1:p) {
    
    x_train <- X_train[, j]
    x_valid <- X_valid[, j]
    
    K_valid_raw <- outer(
      x_valid,
      x_train,
      function(a, b) {
        exp(-gamma * (a - b)^2)
      }
    )
    
    G_valid_list[[j]] <-
      H_valid %*% K_valid_raw %*% H_train
  }
  
  return(G_valid_list)
}


# ============================================================
# Transformed response construction
# ============================================================

# Construct transformed scalar response
compute_Y_vec <- function(Y_list, Y_ref) {
  
  sapply(
    Y_list,
    function(Y_i) {
      cholesky_inner_product(Y_i, Y_ref) -
        cholesky_to_zero(Y_i)
    }
  )
}


# ============================================================
# ADMM Elastic Net solver
# ============================================================

admm_elastic_net <- function(
    G_list,
    Y_vec,
    lambda1,
    lambda2,
    rho = 1,
    max_iter = 1000,
    tol = 1e-6
) {
  
  p <- length(G_list)
  n <- length(Y_vec)
  
  # Initialization
  alpha <- matrix(0, n, p)
  z <- matrix(0, n, p)
  u <- matrix(0, n, p)
  
  # Numerical stabilization
  epsilon <- 1e-4
  
  G_list <- lapply(
    G_list,
    function(G_j) {
      G_j + epsilon * diag(n)
    }
  )
  
  # Precompute inverse matrices
  G_inv_list <- lapply(
    G_list,
    function(G_j) {
      solve(
        (1 / n) * G_j %*% G_j +
          rho * G_j +
          lambda2 * diag(n)
      )
    }
  )
  
  # ADMM iterations
  for (iter in 1:max_iter) {
    
    residual <- Y_vec -
      Reduce(
        `+`,
        lapply(
          1:p,
          function(i) {
            G_list[[i]] %*% alpha[, i]
          }
        )
      )
    
    # Alpha update
    alpha_new <- matrix(0, n, p)
    
    for (j in 1:p) {
      
      rhs <- (1 / n) *
        G_list[[j]] %*%
        (residual + G_list[[j]] %*% alpha[, j]) +
        rho * G_list[[j]] %*% (z[, j] - u[, j])
      
      alpha_new[, j] <- G_inv_list[[j]] %*% rhs
    }
    
    # Group soft-thresholding
    z_new <- matrix(0, n, p)
    
    for (j in 1:p) {
      
      v_j <- alpha_new[, j] + u[, j]
      
      norm_vj <- sqrt(
        sum(v_j * (G_list[[j]] %*% v_j))
      )
      
      z_new[, j] <- max(
        0,
        1 - lambda1[j] / (rho * norm_vj)
      ) * v_j
    }
    
    # Dual update
    u_new <- u + (alpha_new - z_new)
    
    # Convergence check
    if (max(abs(alpha_new - alpha)) < tol) {
      break
    }
    
    alpha <- alpha_new
    z <- z_new
    u <- u_new
  }
  
  return(z)
}


# ============================================================
# Adaptive SCAD refinement
# ============================================================

adaptive_lasso_admm <- function(
    G_list,
    Y_vec,
    lambda1,
    lambda2
) {
  
  # Initial Elastic Net estimator
  alpha_init <- admm_elastic_net(
    G_list,
    Y_vec,
    lambda1 - lambda1,
    lambda2
  )
  
  # Group norms
  value <- sapply(
    1:length(G_list),
    function(j) {
      sum(
        alpha_init[, j] *
          (G_list[[j]] %*% alpha_init[, j])
      )
    }
  )
  
  # Adaptive SCAD weights
  w <- scadrightderv(
    tt = value^0.5,
    lambda = lambda1
  )
  
  # Refit
  return(
    admm_elastic_net(
      G_list,
      Y_vec,
      w,
      lambda2
    )
  )
}


# ============================================================
# Validation-based tuning parameter selection
# ============================================================

select_lambda1 <- function(
    X_valid,
    X_train,
    G_train,
    Y_train,
    G_valid,
    Y_valid,
    lambda1_candidates,
    lambda2,
    method
) {
  
  best_lambda1 <- NULL
  best_mse <- Inf
  
  for (lambda1 in lambda1_candidates) {
    
    # Initial fit
    if (method == "standard") {
      
      result_alpha <- admm_elastic_net(
        G_train,
        Y_train,
        rep(lambda1, length(G_train)),
        lambda2
      )
      
    } else {
      
      result_alpha <- adaptive_lasso_admm(
        G_train,
        Y_train,
        rep(lambda1, length(G_train)),
        lambda2
      )
    }
    
    # Group norms
    value <- sapply(
      1:length(G_train),
      function(j) {
        sum(
          result_alpha[, j] *
            (G_train[[j]] %*% result_alpha[, j])
        )
      }
    )
    
    # Refit selected model
    result_alpha <- admm_elastic_net(
      G_train,
      Y_train,
      (1 - (value > 0)) * 1e8,
      lambda2
    )
    
    # Validation kernel matrices
    G_valid_list <- compute_validation_gram(
      X_valid,
      X_train,
      gamma = 1
    )
    
    # Validation prediction
    predictions <- Reduce(
      `+`,
      lapply(
        1:length(G_train),
        function(j) {
          G_valid_list[[j]] %*%
            result_alpha[, j]
        }
      )
    )
    
    # Validation error
    mse <- mean((Y_valid - predictions)^2)
    
    
    # Update best tuning parameter
    if (Re(mse) < best_mse) {
      
      best_mse <- Re(mse)
      best_lambda1 <- lambda1
      result <- which(value > 0)
    }
  }
  
  return(
    list(
      bl = best_lambda1,
      result = result
    )
  )
}


# ============================================================
# Main simulation experiment
# ============================================================

run_experiment_table <- function(
    n_values,
    p,
    rho,
    lambda1_candidates,
    lambda2,
    num_trials = 100
) {
  
  results <- data.frame()
  
  for (n in n_values) {
    
    print(paste("Running n =", n))
    
    selection_counts_standard <- numeric(p)
    selection_counts_adaptive <- numeric(p)
    
    for (trial in 1:num_trials) {
      
      # ======================================================
      # Training data
      # ======================================================
      
      n_train <- n
      
      X_train <- generate_X(n_train, p, rho)
      
      Y_train_list <- generate_Y_spd_nonlinear(X_train)
      
      X_train <- X_train - mean(X_train)
      
      G_train <- compute_gram_matrices(X_train)
      
      # Reference SPD matrix
      Y_ref_train <- Reduce(
        "+",
        Y_train_list
      ) / length(Y_train_list) + 5
      
      # Construct transformed response
      Y_vec_train <- compute_Y_vec(
        Y_train_list,
        Y_ref_train
      )
      
      sd_y <- sd(Y_vec_train)
      
      Y_vec_train <- (
        Y_vec_train - mean(Y_vec_train)
      ) / sd_y
      
      
      # ======================================================
      # Validation data
      # ======================================================
      
      n_valid <- n
      
      X_valid <- generate_X(n_valid, p, rho)
      
      Y_valid_list <- generate_Y_spd_nonlinear(X_valid)
      
      X_valid <- X_valid - mean(X_valid)
      
      G_valid <- compute_gram_matrices(X_valid)
      
      Y_ref_valid <- Y_ref_train
      
      Y_vec_valid <- compute_Y_vec(
        Y_valid_list,
        Y_ref_valid
      )
      
      Y_vec_valid <- (
        Y_vec_valid - mean(Y_vec_valid)
      ) / sd_y
      
      
      # ======================================================
      # Standard Elastic Net
      # ======================================================
      
      best_lambda1_standard <- select_lambda1(
        X_valid,
        X_train,
        G_train,
        Y_vec_train,
        G_valid,
        Y_vec_valid,
        lambda1_candidates,
        lambda2,
        "standard"
      )
      
      selection_counts_standard[
        best_lambda1_standard$result
      ] <-
        selection_counts_standard[
          best_lambda1_standard$result
        ] + 1
      
      
      # ======================================================
      # Adaptive SCAD refinement
      # ======================================================
      
      best_lambda1_adaptive <- select_lambda1(
        X_valid,
        X_train,
        G_train,
        Y_vec_train,
        G_valid,
        Y_vec_valid,
        lambda1_candidates,
        lambda2,
        "adaptive"
      )
      
      selection_counts_adaptive[
        best_lambda1_adaptive$result
      ] <-
        selection_counts_adaptive[
          best_lambda1_adaptive$result
        ] + 1
    }
    
    # Selection frequencies
    selection_freq_standard <-
      selection_counts_standard / num_trials
    
    selection_freq_adaptive <-
      selection_counts_adaptive / num_trials
    
    
    # Store results
    standard_df <- data.frame(
      n = n,
      Method = "Standard",
      t(selection_freq_standard)
    )
    
    adaptive_df <- data.frame(
      n = n,
      Method = "Adaptive",
      t(selection_freq_adaptive)
    )
    
    results <- rbind(
      results,
      standard_df,
      adaptive_df
    )
  }
  
  colnames(results) <-
    c("n", "Method", paste0("X", 1:p))
  
  print("Final selection frequencies:")
  print(results)
  
  return(results)
}


# ============================================================
# Run simulation
# ============================================================

n_values <- c(200)

results <- run_experiment_table(
  n_values = n_values,
  p = 10,
  rho = 0.5,
  lambda1_candidates = seq(
    0.01,
    0.4,
    length.out = 40
  ),
  lambda2 = 0.1,
  num_trials = 10
)