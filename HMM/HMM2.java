import java.util.Scanner;
import java.util.stream.Stream;

public class HMM2
{
    public static void main(String[] args)
    {
        Scanner in = new Scanner(System.in);

        // Parse data
        // [row][col]
        double[][] transitionMatrix = stringToMatrix(in.nextLine()); // A: Transition probability matrix
        double[][] emissionMatrix = stringToMatrix(in.nextLine()); // B: Observation probability matrix
        double[][] initialProbMatrix = stringToMatrix(in.nextLine()); // pi: Initial probability vector

        int[] emissionSequence = stringToVector(in.nextLine()); // Observation sequence

        in.close();

        // Initialize delta
        int Tr = emissionSequence.length; // T is the number observation sequence elements
        int Nc = transitionMatrix.length; // number of states

        double[][] delta = new double[Tr][Nc];
        int[][] deltaIndices = new int[Tr][Nc];

        //delta

        // first row of delta : delta_1(i) = pi_i * b_i(O_1)
        // here O_1 is the element at index 0 in emissionSequence array
        for (int i = 0; i < Nc; i++) {
            delta[0][i] = initialProbMatrix[0][i] * emissionMatrix[i][emissionSequence[0]];
        }

        // Compute rest of delta : 
        // delta_t(i) = max[ delta_t-1(j) * a_ji * b_i(O_t) ]
        // deltaIDX_t(i) = argmax[ delta_t-1(j) * a_ji * b_i(O_t) ]
        // So, if the algorithm determined delta_t(i) for state i at time t to have been preceded by state j, then deltaIDX_t(i) = j
        // starting from t = 1 because we alreday did first row t = 0
        for (int t = 1; t < Tr; t++) {
            for (int i = 0; i < Nc; i++) {
                double max = -1;
                int max_index = -1;

                for (int j = 0; j < Nc; j++) {
                    double tmp = delta[t - 1][j] * transitionMatrix[j][i] * emissionMatrix[i][emissionSequence[t]];

                    if (tmp > max) {
                        max = tmp;
                        max_index = j;
                    }
                }
                delta[t][i] = max;
                deltaIndices[t][i] = max_index;
            }
        }

        // Compute the most probable sequence of states as zero-based indices
        int[] x_star = new int[Tr];

        // getting the last element of our X* array which is X*_T from the last row of delta matrix
        // X*_T = argmax[ delta_T(j)]
        double max = -1;
        for (int j = 0; j < Nc; j++) {
            double curr = delta[Tr - 1][j];

            if (curr > max) {
                max = curr;
                x_star[Tr - 1] = j; // X*_T = argmax[ delta_T(j)]
            }
        }

        // getting the rest of our X* array by tracing back
        // X*_t = deltaIDX_t+1(X*_t+1)
        // we start from T - 2 because we alredy got the last element of X*
        for (int t = Tr - 2; t >= 0; t--) {
            x_star[t] = deltaIndices[t + 1][x_star[t + 1]];
        }


        

        // Print answer
        StringBuilder sb = new StringBuilder();
        for (int val : x_star) {
            sb.append(val);
            sb.append(' ');
        }
        System.out.println(sb);

    }

    // Helper functions // 

    public static double[][] stringToMatrix(String str)
    {
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
}
