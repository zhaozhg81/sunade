/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sunade.stat.general;

/**
 *
 * @author zhaozhg
 */
public class QUICKSORT {
    public QUICKSORT( ) {    }
   
    public static void quicksort(double[] main, int[] index) {
        quicksort(main, index, 0, index.length - 1);
    }

    // quicksort a[left] to a[right]
    public static void quicksort(double[] a, int[] index, int left, int right) {
        if (right <= left) return;
        int i = partition(a, index, left, right);
        quicksort(a, index, left, i-1);
        quicksort(a, index, i+1, right);
    }

    // partition a[left] to a[right], assumes left < right
    private static int partition(double[] a, int[] index, int left, int right) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (less(a[++i], a[right]))      // find item on left to swap
            ;                               // a[right] acts as sentinel
            while (less(a[right], a[--j]))      // find item on right to swap
                if (j == left) break;           // don't go out-of-bounds
            if (i >= j) break;                  // check if pointers cross
            exch(a, index, i, j);               // swap two elements into place
        }
        exch(a, index, i, right);               // swap with partition element
        return i;
    }

    // is x < y ?
    private static boolean less(double x, double y) {
        return (x < y);
    }

    // exchange a[i] and a[j]
    private static void exch(double[] a, int[] index, int i, int j) {
        double swap = a[i];
        a[i] = a[j];
        a[j] = swap;
        int b = index[i];
        index[i] = index[j];
        index[j] = b;
    }       
    
    
    // quick sort without rearranging the index
    
    public static void quicksort(double[] main) {
        quicksort(main, 0, main.length - 1);
    }

    // quicksort a[left] to a[right]
    public static void quicksort(double[] a, int left, int right) {
        if (right <= left) return;
        int i = partition(a, left, right);
        quicksort(a, left, i-1);
        quicksort(a, i+1, right);
    }

    // partition a[left] to a[right], assumes left < right
    private static int partition(double[] a, int left, int right) {
        int i = left - 1;
        int j = right;
        while (true) {
            while (less(a[++i], a[right]))      // find item on left to swap
            ;                               // a[right] acts as sentinel
            while (less(a[right], a[--j]))      // find item on right to swap
                if (j == left) break;           // don't go out-of-bounds
            if (i >= j) break;                  // check if pointers cross
            exch(a, i, j);               // swap two elements into place
        }
        exch(a, i, right);               // swap with partition element
        return i;
    }


    // exchange a[i] and a[j]
    private static void exch(double[] a, int i, int j) {
        double swap = a[i];
        a[i] = a[j];
        a[j] = swap;
    }       

}
