import java.util.Scanner;
import java.util.stream.Stream;

public class HMM3
{
    public static void main(String[] args)
    {
        Scanner in = new Scanner(System.in);

        // Parse data (1. Initialize ¸ lambda=(A,B,pi))
        // [row][col]
        double[][] transitionMatrix = stringToMatrix(in.nextLine()); // Transition probability matrix
        double[][] emissionMatrix = stringToMatrix(in.nextLine()); // Observation probability matrix
        double[][] initialProbMatrix = stringToMatrix(in.nextLine()); // Initial probability vector

        int[] emissionSequence = stringToVector(in.nextLine()); // Observation sequence

        in.close();

        /*
         *   2. Compute alpha, beta, digamma, gamma values
         *   3. Re-estimate ¸ lambda=(A,B,pi)
         *   4. Repeat from 2. until convergence
         */
        baumWelch(transitionMatrix, emissionMatrix, initialProbMatrix, emissionSequence);
    }

    // The parts related to alpha are similar to HMM1
    public static void alphaPass(double[][] alpha, double[][] transitionMatrix, double[][] emissionMatrix,
                                 double[][] initialProbMatrix, int[] emissionSequence,
                                 double[] c, int Tr, int Nc)
    { 
        // Compute alpha[0][i] : first row of alpha : alpha_1(i) = pi_i * b_i(O_1)
        c[0] = 0; //ST
        for (int i = 0; i < Nc; i++) {
            alpha[0][i] = initialProbMatrix[0][i] * emissionMatrix[i][emissionSequence[0]];
            //ST
            //#### first element of C = the sum of first row of alpha
            c[0] += alpha[0][i]; 
        }

        //####
        //ST 
        // Scale the alpha[0][i]
        c[0] = 1 / c[0];
        for (int i = 0; i < Nc; i++) {
            alpha[0][i] *= c[0]; //alpha_1(i) = alpha_1(i) * c_0
        }

        // Compute a[t][i]
        // Compute rest of alpha : alpha_t(i) = [ sum( alpha_t-1(j) * a_ji ) ] * b_i(O_t)
        //similar to HMM1
        for (int t = 1; t < Tr; t++) {
            c[t] = 0; //ST
            for (int i = 0; i < Nc; i++) {
                alpha[t][i] = 0;
                for (int j = 0; j < Nc; j++) {
                    alpha[t][i] += alpha[t - 1][j] * transitionMatrix[j][i];
                }
                alpha[t][i] *= emissionMatrix[i][emissionSequence[t]];
                //ST
                //#### rest of C = the sum of rest of alpha
                c[t] += alpha[t][i]; 
            }

            //####
            //ST
            // Scale alpha[t][i]
            c[t] = 1 / c[t];
            for (int i = 0; i < Nc; i++) {
                alpha[t][i] *= c[t]; //alpha_t(i) = alpha_t(i) * c_t
            }
        }
    }

    public static void betaPass(double[][] beta, double[][] transitionMatrix, double[][] emissionMatrix,
                                int[] emissionSequence, double[] c,
                                int Tr, int Nc)
    {
        //ST
        // Let beta[Tr - 1][i] = 1, scaled by c[Tr - 1]
        // last row of beta (WITHOUT scalling): beta_T(i) = 1
        for (int i = 0; i < Nc; i++) {
            beta[Tr - 1][i] = c[Tr - 1];
        }

        // Beta-pass
        // Compute rest of beta : (normal) 
        // beta_t(i) = sum[ a_ij * b_j(O_t+1) * beta_t+1(j) ]
        for (int t = Tr - 2; t >= 0; t--) {
            for (int i = 0; i < Nc; i++) {
                beta[t][i] = 0;
                for (int j = 0; j < Nc; j++) {
                    beta[t][i] += transitionMatrix[i][j] * emissionMatrix[j][emissionSequence[t + 1]] * beta[t + 1][j];
                }
                //####
                //ST 
                beta[t][i] *= c[t];  
            }
        }
    }

    public static void diGamma(double[][] alpha, double[][] beta, double[][] gamma, double[][][] diGamma,
                               double[][] transitionMatrix, double[][] emissionMatrix,
                               int[] emissionSequence, int Tr, int Nc)
    {

        // No need to normalize gamma[t][i][j] since using scaled alpha and beta
        // digama_t(i,j) = alpha_t(i) * a_ij * b_j(O_t+1) * beta_t+1(j)
        for (int t = 0; t < Tr - 1; t++) {
            for (int i = 0; i < Nc; i++) {
                gamma[t][i] = 0;
                for (int j = 0; j < Nc; j++) {
                    diGamma[t][i][j] = alpha[t][i] * transitionMatrix[i][j] * emissionMatrix[j][emissionSequence[t + 1]] * beta[t + 1][j];//ST: no need for denominator

                    gamma[t][i] += diGamma[t][i][j]; //ST
                }
            }
        }
        // Special case for gamma[Tr - 1][i] = alpha[Tr - 1][i] (as above, no need to normalize)
        System.arraycopy(alpha[Tr - 1], 0, gamma[Tr - 1], 0, Nc);
    }

    public static void baumWelch(double[][] transitionMatrix, double[][] emissionMatrix,
                                 double[][] initialProbMatrix, int[] emissionSequence)
    {
        //ST
        int iters = 0;
        int maxIters = 100;
        double oldLogProb = Double.NEGATIVE_INFINITY;

        int Tr = emissionSequence.length;
        int Nc = transitionMatrix.length;

        double[] c = new double[emissionSequence.length];
        double[][] alpha = new double[Tr][Nc];
        double[][] beta = new double[Tr][Nc];
        double[][] gamma = new double[Tr][Nc];
        double[][][] diGamma = new double[Tr][Nc][Nc];

        while (true) {
            alphaPass(alpha, transitionMatrix, emissionMatrix, initialProbMatrix, emissionSequence, c, Tr, Nc);
            betaPass(beta, transitionMatrix, emissionMatrix, emissionSequence, c, Tr, Nc);
            diGamma(alpha, beta, gamma, diGamma, transitionMatrix, emissionMatrix, emissionSequence, Tr, Nc);

            // Re-estimate initial probability matrix
            // estimate pi : pi_i = gamma_1(t) -> fisrt row of gamma
            System.arraycopy(gamma[0], 0, initialProbMatrix[0], 0, Nc);


            double denominator;
            double numerator;
            // Re-estimate transition (A) matrix
            // numerator = sum[ digamma_t(i,j)] for t from 1 to T-1
            // denominator = sum[ gamma_t(i)] for t from 1 to T-1
            for (int i = 0; i < Nc; i++) {
                denominator = 0;
                for (int t = 0; t < Tr - 1; t++) {
                    denominator += gamma[t][i];
                }
                for (int j = 0; j < Nc; j++) {
                    numerator = 0;

                    for (int t = 0; t < Tr - 1; t++) {
                        numerator += diGamma[t][i][j];
                    }
                    transitionMatrix[i][j] = numerator / denominator;
                }
            }

            // Re-estimate emission (B) matrix
            // numerator = sum[  1 (if O_t = j) * gamma_t(i) ] for t from 1 to T-1
            // denominator = sum[ gamma_t(i)] for t from 1 to T-1
            for (int i = 0; i < Nc; i++) {
                denominator = 0;
                for (int t = 0; t < Tr; t++) {
                    denominator += gamma[t][i];
                }
                for (int j = 0; j < emissionMatrix[0].length; j++) {
                    numerator = 0;

                    for (int t = 0; t < Tr; t++) {
                        if (emissionSequence[t] == j) {
                            numerator += gamma[t][i];
                        }
                    }
                    emissionMatrix[i][j] = numerator / denominator;
                }
            }

            //ST
            // Compute log[P(O|lambda]] 
            // logProb = logProb + log(c_i)
            double logProb = 0;
            for (int i = 0; i < Tr; i++) {
                logProb += Math.log(c[i]);
            }
            logProb = -logProb;

            //ST
            iters++;
            if (iters < maxIters && logProb > oldLogProb) {
                oldLogProb = logProb;
            } else {
                printAnswer(transitionMatrix);
                printAnswer(emissionMatrix);
                System.exit(0);
            }
        }
    }

    public static double[][] stringToMatrix(String str)
    {
        str = str.trim(); // Remove empty spaces at beginning and end of line

        double[] array = Stream.of(str.split(" ")) // Split line into individual number strings
                .mapToDouble(Double::parseDouble) // Parse strings to doubles
                .toArray(); // Save as array

        // Nr of rows defined as first number in string, column as second
        int rows = (int) array[0];
        int cols = (int) array[1];

        double[][] matrix = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = array[(i * cols) + j + 2]; // +2 to skip row/col definitions at beginning of line
            }
        }
        return matrix;
    }

    public static int[] stringToVector(String str)
    {
        return Stream.of(str.split(" ")) // Split line into individual number strings
                .skip(1) // Skip length definition at beginning of line
                .mapToInt(Integer::parseInt) // Parse strings to integers
                .toArray(); // Return as array
    }

    public static void printAnswer(double[][] ans)
    {
        int rows = ans.length;
        int cols = ans[0].length;

        StringBuilder sb = new StringBuilder();
        sb.append(rows);
        sb.append(' ');
        sb.append(cols);
        sb.append(' ');

        for (double[] row : ans) {
            for (double col : row) {
                sb.append(col);
                sb.append(' ');
            }
        }
        System.out.println(sb);
    }
}