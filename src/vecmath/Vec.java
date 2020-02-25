package vecmath;

import java.io.Serializable;
import java.util.Arrays;
import java.util.function.Function;

public class Vec implements Serializable {
    
    private static final long serialVersionUID = -4924750092823606522L;
    
    final float[] A;
    
    public Vec(float... a) { A = a; }
    
    public Vec(int size, Function<Integer, Float> generator) {
        
        A = new float[size];
        
        for(int i = 0; i < size; i++) {
            A[i] = generator.apply(i);
        }
    }
    
    public Vec add(Vec v) {
        
        if(dimensions() != v.dimensions())
            throw new IllegalArgumentException();
        
        return new Vec(dimensions(), i -> A[i] + v.A[i]);
    }
    
    public Vec sub(Vec v) {
        return add(v.negate());
    }
    
    public Vec mul(float s) {
        return new Vec(dimensions(), i -> A[i] * s);
    }
    
    public Vec div(float s) {
        return mul(1.0F / s);
    }
    
    public Vec negate() {
        return mul(-1.0F);
    }
    
    public float dot(Vec v) {
        
        if(dimensions() != v.dimensions())
            throw new IllegalArgumentException();
        
        float dot = 0.0F;
        
        for(int i = 0; i < dimensions(); i++) {
            dot += A[i] * v.A[i];
        }
        return dot;
    }
    
    public Vec cross(Vec v) {
        
        if(dimensions() != 3 || v.dimensions() != 3)
            throw new IllegalArgumentException();
        
        return new Vec(
            A[1]*v.A[2] - A[2]*v.A[1],
            A[2]*v.A[0] - A[0]*v.A[2],
            A[0]*v.A[1] - A[1]*v.A[0]);
    }
    
    public Mat outer(Vec v) {
        return new Mat(dimensions(), v.dimensions(),
            (r, c) -> A[r] * v.A[c]);
    }
    
    public float triple(Vec v1, Vec v2) {
        return dot(v1.cross(v2));
    }
    
    public float length() {
        return (float) Math.sqrt(lengthSquared());
    }
    
    public float lengthSquared() {
        
        float len = 0.0F;
        
        for(int i = 0; i < dimensions(); i++) {
            len += A[i] * A[i];
        }
        return len;
    }
    
    public Vec setLength(float l) {
        return normalize().mul(l);
    }
    
    public Vec normalize() {
        float len = length();
        if(len == 0) throw new IllegalArgumentException();
        return new Vec(dimensions(), i -> A[i] / len);
    }
    
    public float angle(Vec v) {
        return (float) Math.acos(normalize().dot(v.normalize()));
    }
    
    public float distance(Vec v) {
        return sub(v).length();
    }
    
    public Vec project(Vec v) {
        return v.normalize().mul(dot(v.normalize()));
    }
    
    public Vec reflect(Vec v) {
        Vec proj = project(v);
        return sub(proj).negate().add(proj);
    }
    
    public Vec lerp(Vec v, float s) {
        return mul(1.0F-s).add(v.mul(s));
    }
    
    public Vec midpoint(Vec v) {
        return lerp(v, 0.5F);
    }
    
    public Vec append(Vec v) {
        return new Vec(dimensions() + v.dimensions(),
            i -> i<dimensions() ? A[i] : v.A[i-dimensions()]);
    }
    
    public Vec toHomogenous() {
        return append(new Vec(1.0F));
    }
    
    public Vec toCartesian() {
        return new Vec(dimensions()-1,
            i -> A[i] / A[dimensions()-1]);
    }
    
    public Vec changeBasis(Mat basis) {
        return basis.invert().mul(this);
    }
    
    public float val(int i) {
        return A[i];
    }
    
    public int dimensions() {
        return A.length;
    }
    
    @Override
    public String toString() {
        return Arrays.toString(A);
    }
    
    public static Vec zero(int dim) {
        return repeat(0.0F, dim);
    }
    
    public static Vec one(int dim) {
        return repeat(1.0F, dim);
    }
    
    public static Vec repeat(float val, int dim) {
        return new Vec(dim, i -> val);
    }
    
    public static Vec basis(int axis, int dim) {
        return new Vec(dim, i -> i==axis ? 1.0F : 0.0F);
    }
}