package vecmath;

import java.io.Serializable;
import java.util.function.Function;

public class Quat implements Serializable {
    
    private static final long serialVersionUID = 7322359333184172345L;
    
    final float[] A;
    
    public Quat(float[] a) { A = a; }
    
    public Quat(Vec v) {
        this(0, v);
    }
    
    public Quat(float s, Vec v) {
        this(i -> i==0 ? s : v.A[i-1]);
        if(v.dimensions() != 3) throw new IllegalArgumentException();
    }
    
    public Quat(Function<Integer, Float> generator) {
        
        A = new float[4];
        
        for(int i = 0; i < 4; i++) {
            A[i] = generator.apply(i);
        }
    }
    
    public Quat add(Quat q) {
        return new Quat(i -> A[i] + q.A[i]);
    }
    
    public Quat sub(Quat q) {
        return add(q.negate());
    }
    
    public Quat mul(Quat q) {
        
        float s1 = scalar(), s2 = q.scalar();
        Vec v1 = vector(), v2 = q.vector();
        
        return new Quat(s1*s2 + v1.dot(v2),
            v2.mul(s1).add(v1.mul(s2)).add(v1.cross(v2)));
    }
    
    public Vec mul(Vec v) {
        return mul(new Quat(v)).vector();
    }
    
    public Quat mul(float s) {
        return new Quat(i -> A[i] * s);
    }
    
    public Quat div(float s) {
        return mul(1.0F / s);
    }
    
    public Quat negate() {
        return mul(-1.0F);
    }
    
    public float dot(Quat v) {
        
        float dot = 0.0F;
        
        for(int i = 0; i < 4; i++) {
            dot += A[i] * v.A[i];
        }
        return dot;
    }
    
    public Quat conjugate() {
        return new Quat(i -> i==0 ? A[i] : -A[i]);
    }
    
    public Quat invert() {
        return conjugate().div(norm());
    }
    
    public float norm() {
        return (float) Math.sqrt(scalar()*scalar()
            + vector().lengthSquared());
    }
    
    public Quat normalize() {
        float norm = norm();
        if(norm == 0) throw new IllegalArgumentException();
        return div(norm);
    }
    
    public Vec rotate(Vec v) {
        return mul(new Quat(v)).mul(invert()).vector();
    }
    
    public Vec axis() {
        return vector().normalize();
    }
    
    public float angle() {
        return 2.0F * (float) Math.atan2(vector().length(), scalar());
    }
    
    public Quat slerp(Quat q, float t) {
        
        q = q.normalize();
        float dot = dot(q);
        if(dot < 0.0F) {
            q = q.negate();
            dot *= -1.0F;
        }
        
        if(dot > 0.9995F) {
            return mul(1.0F-t).add(q.mul(t)).normalize();
        }
        
        float a = (float) Math.acos(dot);
        float s1 = (float) (Math.sin(a*t) / Math.sin(a));
        float s0 = (float) Math.cos(a*t) - dot * s1;
        
        return mul(s0).add(q.mul(s1)).normalize();
    }
    
    public float scalar() {
        return A[0];
    }
    
    public Vec vector() {
        return new Vec(3, i -> A[i+1]);
    }
    
    public float val(int i) {
        return A[i];
    }
    
    public static Quat identity() {
        return new Quat(1.0F, Vec.zero(3));
    }
    
    public static Quat angleAxis(float angle, Vec axis) {
        return new Quat((float) Math.cos(angle / 2.0F),
            axis.setLength((float) Math.sin(angle / 2.0F)));
    }
    
    public static Quat angleAxis(float angle, int axis) {
        return angleAxis(angle, Vec.basis(axis, 3));
    }
    
    public static Quat euler(float yaw, float pitch, float roll) {
        return angleAxis(yaw, Vec.basis(1, 3))
            .mul(angleAxis(-pitch, Vec.basis(0, 3)))
            .mul(angleAxis(roll, Vec.basis(2, 3)));
    }
}