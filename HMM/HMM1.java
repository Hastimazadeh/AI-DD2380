import java.util.Scanner;
import java.util.stream.Stream;

public class HMM1
{
    public static void main(String[] args)
    {
        Scanner in = new Scanner(System.in);

        // Parse data
        // [row][col]
        double[][] transitionMatrix = stringToMatrix(in.nextLine()); // A : Transition probability matrix
        double[][] emissionMatrix = stringToMatrix(in.nextLine()); // B: Observation probability matrix
        double[][] initialProbMatrix = stringToMatrix(in.nextLine()); // pi: Initial probability vector

        int[] emissionSequence = stringToVector(in.nextLine()); // Observation sequence

        in.close();

        // Initialize alpha
        int Tr = emissionSequence.length; // T is the number observation sequence elements
        int Nc = transitionMatrix.length; // number of states

        double[][] alpha = new double[Tr][Nc];

        // alpha 

        // first row of alpha : alpha_1(i) = pi_i * b_i(O_1)
        // here O_1 is the element at index 0 in emissionSequence array
        for (int i = 0; i < Nc; i++) {
            alpha[0][i] = initialProbMatrix[0][i] * emissionMatrix[i][emissionSequence[0]];
        }

        // Compute rest of alpha : alpha_t(i) = [ sum( alpha_t-1(j) * a_ji ) ] * b_i(O_t)
        // starting from t = 1 because we alreday did first row t = 0
        for (int t = 1; t < Tr; t++) {
            for (int i = 0; i < Nc; i++) {
                double sum = 0;

                for (int j = 0; j < Nc; j++) {
                    sum += alpha[t - 1][j] * transitionMatrix[j][i];
                }
                alpha[t][i] = sum * emissionMatrix[i][emissionSequence[t]];
            }
        }


        // Compute the probability of having observed the given observation sequence O1:T 
        // = sum( alpha_T(i)) where i is 1 - N
        // alpha_T is the last row of our alpha matrix
        // Compute the probability of the given sequence as a single scalar
        double answer = 0;
        for (int i = 0; i < Nc; i++) {
            answer += alpha[Tr - 1][i];
        }
        System.out.println(answer);
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
