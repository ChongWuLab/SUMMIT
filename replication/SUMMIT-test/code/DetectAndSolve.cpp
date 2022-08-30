#include <armadillo>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

int MCPDetect(const arma::vec &r, const arma::mat &LD, arma::vec lambda, const double s, const double gamma, const int max_iteration, const double threshold, int index_lambda)
{
    const int n_var = r.size();
    double lambda_used, z, dlx;

    arma::mat temp(1, 1);

    arma::vec beta;
    beta.zeros(n_var);

    arma::vec beta_current;
    beta_current.zeros(n_var);

    int problem = 1;

    index_lambda = index_lambda - 1;
    lambda_used = lambda(index_lambda);

    for (int j = 0; j < max_iteration; j++)
    {
        for (int q = 0; q < n_var; q++)
        {
            temp = r(q) - LD.col(q).t() * beta_current + LD(q, q) * beta_current(q);
            z = temp(0, 0);
            if (beta_current(q) <= gamma * lambda_used)
            {
                if (z > lambda_used)
                {
                    beta_current(q) = (z - lambda_used) / (1 + s - 1 / gamma);
                }
                else if (z < -lambda_used)
                {
                    beta_current(q) = (z + lambda_used) / (1 + s - 1 / gamma);
                }
                else
                {
                    beta_current(q) = 0.0;
                }
            }
            else
            {
                beta_current(q) = z / (1 + s);
            }
        }

        dlx = arma::norm(beta_current, 2);

        if (dlx > threshold)
        {
            problem = 2;
            break;
        }

        beta = beta_current;
    }

    return (problem);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat MCP(const arma::vec &r, const arma::mat &LD, arma::vec lambda, const double s, const double gamma, const int max_iteration, const double threshold)
{
    const int n_var = r.size();
    const int n_lambda = lambda.size();
    double lambda_used, z, dlx, dlx2;

    arma::mat res(n_var, n_lambda);
    res.fill(0);

    arma::mat temp(1, 1);

    arma::vec beta;
    beta.zeros(n_var);

    arma::vec beta_current;
    beta_current.zeros(n_var);

    int problem = 0;

    for (int i = 0; i < n_lambda; i++)
    {
        beta_current = beta;
        lambda_used = lambda(i);
        for (int j = 0; j < max_iteration; j++)
        {
            for (int q = 0; q < n_var; q++)
            {
                temp = r(q) - LD.col(q).t() * beta_current + LD(q, q) * beta_current(q);
                z = temp(0, 0);
                if (beta_current(q) <= gamma * lambda_used)
                {
                    if (z > lambda_used)
                    {
                        beta_current(q) = (z - lambda_used) / (1 + s - 1 / gamma);
                    }
                    else if (z < -lambda_used)
                    {
                        beta_current(q) = (z + lambda_used) / (1 + s - 1 / gamma);
                    }
                    else
                    {
                        beta_current(q) = 0.0;
                    }
                }
                else
                {
                    beta_current(q) = z / (1 + s);
                }
            }

            arma::vec betaDiff = abs(beta_current - beta);
            dlx = betaDiff.max();
            dlx2 = arma::norm(beta_current, 2);

            if (dlx < threshold)
            {
                break;
            }

            if (dlx2 > 100.0)
            {
                problem = 1;
                break;
            }

            beta = beta_current;
        }

        if (problem == 1)
        {
            break;
        }

        res.col(i) = beta_current;
    }

    return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int ElNetDetect(const arma::vec &r, const arma::mat &LD, const arma::vec lambda, const double s, const double alpha, const int max_iteration, const double threshold, int index_lambda)
{
    const int n_var = r.size();
    double z, lambda_used, dlx;

    arma::mat temp(1, 1);

    arma::vec beta;
    beta.zeros(n_var);

    arma::vec beta_current;
    beta_current.zeros(n_var);

    int problem = 1;

    beta_current = beta;
    index_lambda = index_lambda - 1;
    lambda_used = lambda(index_lambda);

    for (int j = 0; j < max_iteration; j++)
    {
        for (int q = 0; q < n_var; q++)
        {
            temp = r(q) - LD.col(q).t() * beta_current + LD(q, q) * beta_current(q);
            z = temp(0, 0);
            if (z > alpha * lambda_used)
            {
                beta_current(q) = (z - alpha * lambda_used) / (1 + s + 2 * lambda_used * (1 - alpha));
            }
            else if (z < -alpha * lambda_used)
            {
                beta_current(q) = (z + alpha * lambda_used) / (1 + s + 2 * lambda_used * (1 - alpha));
            }
            else
            {
                beta_current(q) = 0.0;
            }
        }

        dlx = arma::norm(beta_current, 2);

        if (dlx > threshold)
        {
            problem = 2;
            break;
        }

        beta = beta_current;
    }

    return (problem);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat ElNet(const arma::vec &r, const arma::mat &LD, const arma::vec lambda, const double s, const double alpha, const int max_iteration, const double threshold)
{
    const int n_var = r.size();
    const int n_lambda = lambda.size();
    double z, lambda_used, dlx, dlx2;

    arma::mat res(n_var, n_lambda);
    res.fill(0);

    arma::mat temp(1, 1);

    arma::vec beta;
    beta.zeros(n_var);

    arma::vec beta_current;
    beta_current.zeros(n_var);

    int problem = 0;

    for (int i = 0; i < n_lambda; i++)
    {
        beta_current = beta;
        lambda_used = lambda(i);

        for (int j = 0; j < max_iteration; j++)
        {
            for (int q = 0; q < n_var; q++)
            {
                temp = r(q) - LD.col(q).t() * beta_current + LD(q, q) * beta_current(q);
                z = temp(0, 0);
                if (z > alpha * lambda_used)
                {
                    beta_current(q) = (z - alpha * lambda_used) / (1 + s + 2 * lambda_used * (1 - alpha));
                }
                else if (z < -alpha * lambda_used)
                {
                    beta_current(q) = (z + alpha * lambda_used) / (1 + s + 2 * lambda_used * (1 - alpha));
                }
                else
                {
                    beta_current(q) = 0.0;
                }
            }

            arma::vec betaDiff = abs(beta_current - beta);
            dlx = betaDiff.max();
            dlx2 = arma::norm(beta_current, 2);

            if (dlx < threshold)
            {
                break;
            }

            if (dlx2 > 100.0)
            {
                problem = 1;
                break;
            }

            beta = beta_current;
        }

        if (problem == 1)
        {
            break;
        }

        res.col(i) = beta_current;
    }

    return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int MNetDetect(const arma::vec &r, const arma::mat &LD, arma::vec lambda, const double s, const double alpha, const double gamma, const int max_iteration, const double threshold, int index_lambda)
{
    const int n_var = r.size();
    double z, lambda_used, dlx, lambda_1, lambda_2;

    arma::mat temp(1, 1);

    arma::vec beta;
    beta.zeros(n_var);

    arma::vec beta_current;
    beta_current.zeros(n_var);

    int problem = 1;

    beta_current = beta;
    index_lambda = index_lambda - 1;
    lambda_used = lambda(index_lambda);
    lambda_1 = lambda_used * alpha;
    lambda_2 = lambda_used * (1 - alpha);

    for (int j = 0; j < max_iteration; j++)
    {
        for (int q = 0; q < n_var; q++)
        {
            temp = r(q) - LD.col(q).t() * beta_current + LD(q, q) * beta_current(q);
            z = temp(0, 0);
            if (abs(z) <= gamma * lambda_1 * (1 + lambda_2))
            {
                if (z > lambda_1)
                {
                    beta_current(q) = (z - lambda_1) / (1 + s + lambda_2 - 1 / gamma);
                }
                else if (z < -lambda_1)
                {
                    beta_current(q) = (z + lambda_1) / (1 + s + lambda_2 - 1 / gamma);
                }
                else
                {
                    beta_current(q) = 0.0;
                }
            }
            else
            {
                beta_current(q) = z / (1 + s + lambda_2);
            }
        }

        dlx = arma::norm(beta_current, 2);

        if (dlx > threshold)
        {
            problem = 2;
            break;
        }

        beta = beta_current;
    }

    return (problem);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat MNet(const arma::vec &r, const arma::mat &LD, arma::vec lambda, const double s, const double alpha, const double gamma, const int max_iteration, const double threshold)
{
    const int n_var = r.size();
    const int n_lambda = lambda.size();
    double z, lambda_used, dlx, dlx2, lambda_1, lambda_2;

    arma::mat res(n_var, n_lambda);
    res.fill(0);

    arma::mat temp(1, 1);

    arma::vec beta;
    beta.zeros(n_var);

    arma::vec beta_current;
    beta_current.zeros(n_var);

    int problem = 0;

    for (int i = 0; i < n_lambda; i++)
    {
        beta_current = beta;
        lambda_used = lambda(i);
        lambda_1 = lambda_used * alpha;
        lambda_2 = lambda_used * (1 - alpha);
        for (int j = 0; j < max_iteration; j++)
        {

            for (int q = 0; q < n_var; q++)
            {
                temp = r(q) - LD.col(q).t() * beta_current + LD(q, q) * beta_current(q);
                z = temp(0, 0);
                if (abs(z) <= gamma * lambda_1 * (1 + lambda_2))
                {
                    if (z > lambda_1)
                    {
                        beta_current(q) = (z - lambda_1) / (1 + s + lambda_2 - 1 / gamma);
                    }
                    else if (z < -lambda_1)
                    {
                        beta_current(q) = (z + lambda_1) / (1 + s + lambda_2 - 1 / gamma);
                    }
                    else
                    {
                        beta_current(q) = 0.0;
                    }
                }
                else
                {
                    beta_current(q) = z / (1 + s + lambda_2);
                }
            }

            arma::vec betaDiff = abs(beta_current - beta);
            dlx = betaDiff.max();

            dlx2 = arma::norm(beta_current, 2);

            if (dlx < threshold)
            {
                break;
            }

            if (dlx2 > 100.0)
            {
                problem = 1;
                break;
            }

            beta = beta_current;
        }

        if (problem == 1)
        {
            break;
        }

        res.col(i) = beta_current;
    }

    return (res);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

int SCADDetect(const arma::vec &r, const arma::mat &LD, arma::vec lambda, const double s, const double gamma, const int max_iteration, const double threshold, int index_lambda)
{
    const int n_var = r.size();
    double lambda_used, z, dlx;

    arma::mat temp(1, 1);

    arma::vec beta;
    beta.zeros(n_var);

    arma::vec beta_current;
    beta_current.zeros(n_var);

    int problem = 1;

    beta_current = beta;
    index_lambda = index_lambda - 1;
    lambda_used = lambda(index_lambda);

    for (int j = 0; j < max_iteration; j++)
    {
        for (int q = 0; q < n_var; q++)
        {
            temp = r(q) - LD.col(q).t() * beta_current + LD(q, q) * beta_current(q);
            z = temp(0, 0);
            if (beta_current(q) <= gamma * lambda_used)
            {
                if (beta_current(q) <= lambda_used)
                {
                    if (z > lambda_used)
                    {
                        beta_current(q) = (z - lambda_used) / (1 + s);
                    }
                    else if (z < -lambda_used)
                    {
                        beta_current(q) = (z + lambda_used) / (1 + s);
                    }
                    else
                    {
                        beta_current(q) = 0.0;
                    }
                }
                else
                {
                    if (z > gamma / (gamma - 1) * lambda_used)
                    {
                        beta_current(q) = (z - gamma / (gamma - 1) * lambda_used) / (1 + s - 1 / (gamma - 1));
                    }
                    else if (z < -gamma / (gamma - 1) * lambda_used)
                    {
                        beta_current(q) = (z + gamma / (gamma - 1) * lambda_used) / (1 + s - 1 / (gamma - 1));
                    }
                    else
                    {
                        beta_current(q) = 0.0;
                    }
                }
            }
            else
            {
                beta_current(q) = z / (1 + s);
            }
        }

        dlx = arma::norm(beta_current, 2);

        if (dlx > threshold)
        {
            problem = 2;
            break;
        }

        beta = beta_current;
    }

    return (problem);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat SCAD(const arma::vec &r, const arma::mat &LD, arma::vec lambda, const double s, const double gamma, const int max_iteration, const double threshold)
{
    const int n_var = r.size();
    const int n_lambda = lambda.size();
    double lambda_used, z, dlx, dlx2;

    arma::mat res(n_var, n_lambda);
    res.fill(0);

    arma::mat temp(1, 1);

    arma::vec beta;
    beta.zeros(n_var);

    arma::vec beta_current;
    beta_current.zeros(n_var);

    int problem = 0;

    for (int i = 0; i < n_lambda; i++)
    {
        beta_current = beta;
        lambda_used = lambda(i);
        for (int j = 0; j < max_iteration; j++)
        {
            for (int q = 0; q < n_var; q++)
            {
                temp = r(q) - LD.col(q).t() * beta_current + LD(q, q) * beta_current(q);
                z = temp(0, 0);
                if (beta_current(q) <= gamma * lambda_used)
                {
                    if (beta_current(q) <= lambda_used)
                    {
                        if (z > lambda_used)
                        {
                            beta_current(q) = (z - lambda_used) / (1 + s);
                        }
                        else if (z < -lambda_used)
                        {
                            beta_current(q) = (z + lambda_used) / (1 + s);
                        }
                        else
                        {
                            beta_current(q) = 0.0;
                        }
                    }
                    else
                    {
                        if (z > gamma / (gamma - 1) * lambda_used)
                        {
                            beta_current(q) = (z - gamma / (gamma - 1) * lambda_used) / (1 + s - 1 / (gamma - 1));
                        }
                        else if (z < -gamma / (gamma - 1) * lambda_used)
                        {
                            beta_current(q) = (z + gamma / (gamma - 1) * lambda_used) / (1 + s - 1 / (gamma - 1));
                        }
                        else
                        {
                            beta_current(q) = 0.0;
                        }
                    }
                }
                else
                {
                    beta_current(q) = z / (1 + s);
                }
            }

            arma::vec betaDiff = abs(beta_current - beta);
            dlx = betaDiff.max();
            dlx2 = arma::norm(beta_current, 2);

            if (dlx < threshold)
            {
                break;
            }

            if (dlx2 > 100.0)
            {
                problem = 1;
                break;
            }

            beta = beta_current;
        }

        if (problem == 1)
        {
            break;
        }

        res.col(i) = beta_current;
    }

    return (res);
}