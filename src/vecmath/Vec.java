package vecmath;

import java.io.Serializable;
import java.util.Arrays;
import java.util.function.Function;

/**
 * Arbitrary-dimensional immutable vector type for computer graphics.
 * Used to represent positions and directions.
 * @author Alec Dorrington
 */
public class Vec implements Serializable {
    
    private static final long serialVersionUID = -4924750092823606522L;
    
    /** The contents of the vector. */
    final float[] A;
    
    /**
     * Construct a new a-dimensional vector.
     * @param a the contents of the vector.
     */
    public Vec(float... a) {
        this(a.length, i -> a[i]);
    }
    
    /**
     * Construct a new dim-dimensional vector using a generator.
     * The generator accepts an index and should return a corresponding value.
     * @param dim the number of dimensions of the vector.
     * @param generator the function for generating the vector's contents.
     */
    public Vec(int dim, Function<Integer, Float> generator) {
        
        A = new float[dim];
        
        for(int i = 0; i < dim; i++) {
            A[i] = generator.apply(i);
        }
    }
    
    /**
     * Get the vector sum of this vector and another.
     * @param v the vector to add.
     * @return the sum of this vector and v.
     */
    public Vec add(Vec v) {
        
        if(dimensions() != v.dimensions())
            throw new IllegalArgumentException();
        
        return new Vec(dimensions(), i -> A[i] + v.A[i]);
    }
    
    /**
     * Get the vector difference of this vector and another.
     * @param v the vector to subtract.
     * @return v subtracted from this vector.
     */
    public Vec sub(Vec v) {
        return add(v.negate());
    }
    
    /**
     * Get the product of this vector and a scalar.
     * @param s the scale factor.
     * @return this vector scaled by s.
     */
    public Vec mul(float s) {
        return new Vec(dimensions(), i -> A[i] * s);
    }
    
    /**
     * Get the quotient of this vector and a scalar.
     * @param s the divisor.
     * @return this vector scaled by 1/s.
     */
    public Vec div(float s) {
        return mul(1.0F / s);
    }
    
    /**
     * Get this vector, negated.
     * @return this vector scaled by -1.
     */
    public Vec negate() {
        return mul(-1.0F);
    }
    
    /**
     * Get the dot product/inner product of this vector and another.
     * An indication as to whether the vectors point in similar directions.
     * @param v the vector to take the dot product with.
     * @return the dot product of this vector and v.
     */
    public float dot(Vec v) {
        
        if(dimensions() != v.dimensions())
            throw new IllegalArgumentException();
        
        float dot = 0.0F;
        
        for(int i = 0; i < dimensions(); i++) {
            dot += A[i] * v.A[i];
        }
        return dot;
    }
    
    /**
     * Get the cross product of this vector and another.
     * The result is a new vector perpendicular to both other vectors.
     * Requires that both vectors are 3-dimensional.
     * @param v the (right) vector to take the cross product with.
     * @return the cross product of this vector and v.
     */
    public Vec cross(Vec v) {
        
        if(dimensions() != 3 || v.dimensions() != 3)
            throw new IllegalArgumentException();
        
        return new Vec(
            A[1]*v.A[2] - A[2]*v.A[1],
            A[2]*v.A[0] - A[0]*v.A[2],
            A[0]*v.A[1] - A[1]*v.A[0]);
    }
    
    /**
     * Get the outer product of this vector and another.
     * Returns the matrix which results from matrix-multiplying the
     * first column vector with the second row vector.
     * @param v the (right) vector to take the outer product with.
     * @return the outer product of this vector and v.
     */
    public Mat outer(Vec v) {
        return new Mat(dimensions(), v.dimensions(),
            (r, c) -> A[r] * v.A[c]);
    }
    
    /**
     * Get the vector triple product of 3 vectors.
     * Defined by taking the cross product of 2 vectors,
     * and then finding the dot product of the result with the third vector.
     * @param v1 the second vector.
     * @param v2 the third vector.
     * @return the triple product of this vector, v1 and v2.
     */
    public float triple(Vec v1, Vec v2) {
        return dot(v1.cross(v2));
    }
    
    /**
     * Get the Euclidean norm (length) of this vector.
     * @return the length of this vector.
     */
    public float length() {
        return (float) Math.sqrt(lengthSquared());
    }
    
    /**
     * Get the square of the Eucliean norm (length) of this vector.
     * Used in cases where cardinal length is not required,
     * as this function is more efficient than calling length().
     * @return the length of this vector, squared.
     */
    public float lengthSquared() {
        
        float len = 0.0F;
        
        for(int i = 0; i < dimensions(); i++) {
            len += A[i] * A[i];
        }
        return len;
    }
    
    /**
     * Get a vector like this but with a particular length.
     * @param len the Eucliean norm (length) of the new vector.
     * @return a new vector of this length.
     */
    public Vec setLength(float len) {
        return normalize().mul(len);
    }
    
    /**
     * Get a unit vector in the same direction as this vector.
     * @return new vector of length 1.
     */
    public Vec normalize() {
        float len = length();
        if(len == 0) throw new IllegalArgumentException();
        return new Vec(dimensions(), i -> A[i] / len);
    }
    
    /**
     * Get the angle between this vector and another.
     * @param v the other vector to measure the angle between.
     * @return the angle, in radians.
     */
    public float angle(Vec v) {
        return (float) Math.acos(normalize().dot(v.normalize()));
    }
    
    /**
     * Get the distance between this vector and another.
     * @param v the other vector to measure the distance between.
     * @return the distance.
     */
    public float distance(Vec v) {
        return sub(v).length();
    }
    
    /**
     * Project this vector onto another.
     * That is, find the component of this vector in the direction of v.
     * @param v the vector to project onto,
     * @return the vector projection.
     */
    public Vec project(Vec v) {
        return v.normalize().mul(dot(v.normalize()));
    }
    
    /**
     * Reflect this vector over the given normal
     * by negating the component parallel to the normal.
     * @param v the normal vector of the surface to reflect off.
     * @return the reflected vector.
     */
    public Vec reflect(Vec v) {
        Vec proj = project(v);
        return sub(proj).negate().add(proj);
    }
    
    /**
     * Linearly-interpolate between this vector and another.
     * Get a new vector somewhere between this vector and v.
     * @param v the vector to interpolate to.
     * @param t the interpolation progress (0.0-1.0).
     * @return the interpolated vector.
     */
    public Vec lerp(Vec v, float t) {
        return mul(1.0F-t).add(v.mul(t));
    }
    
    /**
     * Find the midpoint of this vector and another.
     * @param v the other vector to find the midpoint between.
     * @return the midpoint of the 2 vectors.
     */
    public Vec midpoint(Vec v) {
        return lerp(v, 0.5F);
    }
    
    /**
     * Append another vector to this one, increasing the number of dimensions.
     * @param v the vector to append.
     * @return a new vector containing all the components of both other vectors.
     */
    public Vec append(Vec v) {
        return new Vec(dimensions() + v.dimensions(),
            i -> i<dimensions() ? A[i] : v.A[i-dimensions()]);
    }
    
    /**
     * Convert this vector from cartesian to homogenous coordinates.
     * Increases the dimensionaliy by 1 by adding a divisor term.
     * Allows for the use of general projective transformations.
     * Allows for the representation of points at infinity.
     * @return this vector represented in homogenous coordinates.
     */
    public Vec toHomogenous() {
        return append(new Vec(1.0F));
    }
    
    /**
     * Convert this vector from homogenous to cartesian coordinates.
     * Decreases the dimensionality by 1 by dividing by the last term.
     * Allows for the use of standard vector arithmetic.
     * @return this vector represented in cartesian coordinates.
     */
    public Vec toCartesian() {
        return new Vec(dimensions()-1,
            i -> A[i] / A[dimensions()-1]);
    }
    
    /**
     * Apply a change of basis matrix to this vector,
     * so that it can be represented in a new coordinate system.
     * @param basis a change of basis matrix.
     * @return this vector in the new coordinate system.
     */
    public Vec changeBasis(Mat basis) {
        return basis.invert().mul(this);
    }
    
    /**
     * Get the value of the dim-th component of this vector.
     * @param dim the index of the component to get.
     * @return the value of this component.
     */
    public float val(int dim) {
        return A[dim];
    }
    
    /**
     * Get the number of dimensions present in this vector.
     * @return the dimensionality of this vector.
     */
    public int dimensions() {
        return A.length;
    }
    
    @Override
    public String toString() {
        return Arrays.toString(A);
    }
    
    /**
     * Get a vector containing only zeroes.
     * @param dim the dimensionality of the vector.
     * @return a dim-dimensional vector containing only zeroes.
     */
    public static Vec zero(int dim) {
        return repeat(0.0F, dim);
    }
    
    /**
     * Get a vector containing only ones.
     * @param dim the dimensionality of the vector.
     * @return a dim-dimensional vector containing only ones.
     */
    public static Vec one(int dim) {
        return repeat(1.0F, dim);
    }
    
    /**
     * Get a vector consisting of the same value repeated in each component.
     * @param val the value of each component in the vector.
     * @param dim the dimensionality of the vector.
     * @return a dim-dimensional vector containing only val.
     */
    public static Vec repeat(float val, int dim) {
        return new Vec(dim, i -> val);
    }
    
    /**
     * Return the axis-th basis vector in dim-dimensional space.
     * That is, a vector with 1 in the axis-th position, and 0 elsewhere.
     * @param axis the index of the basis vector to return.
     * @param dim the dimensionality of the space.
     * @return a basis vector.
     */
    public static Vec basis(int axis, int dim) {
        return new Vec(dim, i -> i==axis ? 1.0F : 0.0F);
    }
}