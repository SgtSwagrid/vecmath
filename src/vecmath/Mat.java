package vecmath;

import java.io.Serializable;
import java.util.Arrays;
import java.util.function.BiFunction;

public class Mat implements Serializable {
    
    private static final long serialVersionUID = 3888722482803464661L;
    
    final float[][] A;
    
    public Mat(float[][] a) { A = a; }
    
    public Mat(int height, int width,
            BiFunction<Integer, Integer, Float> generator) {
        
        A = new float[height][width];
        
        for(int r = 0; r < height; r++) {
            for(int c = 0; c < width; c++) {
                A[r][c] = generator.apply(r, c);
            }
        }
    }
    
    public Mat add(Mat m) {
        
        if(height() != m.height() || width() != m.width())
            throw new IllegalArgumentException();
        
        return new Mat(height(), width(),
            (r, c) -> A[r][c] + m.A[r][c]);
    }
    
    public Mat sub(Mat m) {
        return add(m.negate());
    }
    
    public Mat mul(Mat m) {
        
        if(width() != m.height())
            throw new IllegalArgumentException();
        
        return new Mat(height(), m.width(),
            (r, c) -> row(r).dot(m.col(c)));
    }
    
    public Vec mul(Vec v) {
        
        if(width() == v.dimensions()) {
            return new Vec(height(), i -> v.dot(row(i)));
            
        } else if(width() == v.dimensions()+1) {
            Vec vh = v.toHomogenous();
            return new Vec(height(), i -> vh.dot(row(i)))
                .toCartesian();
            
        } else throw new IllegalArgumentException();
    }
    
    public Mat mul(float s) {
        return new Mat(height(), width(),
            (r, c) -> A[r][c] * s);
    }
    
    public Mat div(float s) {
        return mul(1.0F / s);
    }
    
    public Mat negate() {
        return mul(-1.0F);
    }
    
    public Mat power(int i) {
        
        if(i == 1) return this;
        else if(i == 0) return identity(height());
        else if(i < 0) return invert().power(-i);
        else return mul(power(i-1));
    }
    
    public Mat transpose() {
        return new Mat(width(), height(),
            (r, c) -> A[c][r]);
    }
    
    public float determinant() {
        
        if(height() != width())
            throw new IllegalArgumentException();
        
        if(height()==1 && width()==1) {
            return A[0][0];
            
        } else {
            float det = 0.0F;
            for(int c = 0; c < width(); c++) {
                det += A[0][c] * cofactor(0, c);
            }
            return det;
        }
    }
    
    public Mat invert() {
        return adjugate().div(determinant());
    }
    
    public Mat submatrix(int r, int c) {
        return new Mat(height()-1, width()-1,
            (rr, cc) -> A[rr<r ? rr : rr+1][cc<c ? cc : cc+1]);
    }
    
    public float minor(int r, int c) {
        return submatrix(r, c).determinant();
    }
    
    public Mat minor() {
        return new Mat(height(), width(),
            (r, c) -> minor(r, c));
    }
    
    public float cofactor(int r, int c) {
        return minor(r, c) * ((r+c)%2==0 ? 1.0F : -1.0F);
    }
    
    public Mat cofactor() {
        return new Mat(height(), width(),
            (r, c) -> cofactor(r, c));
    }
    
    public Mat adjugate() {
        return cofactor().transpose();
    }
    
    public Mat changeBasis(Mat basis) {
        return basis.invert().mul(this);
    }
    
    public Mat globalTransform(Mat m) {
        return m.mul(this);
    }
    
    public Mat extrinsicTransform(Mat m) {
        Mat tm = translate(getTranslation());
        return tm.mul(m).mul(tm.invert()).mul(this);
    }
    
    public Mat intrinsicTransform(Mat m) {
        return translate(getTranslation()).mul(getRotation())
            .mul(m).mul(scale(getScale()));
    }
    
    public Vec getTranslation() {
        return col(width()-1).toCartesian();
    }
    
    public Mat getRotation() {
        return globalTransform(scale(getScale()).invert()
            .mul(translate(getTranslation().negate())));
    }
    
    public Vec getScale() {
        return new Vec(width()-1, i -> col(i).length());
    }
    
    public Vec row(int r) {
        return new Vec(A[r]);
    }
    
    public Vec col(int c) {
        return new Vec(height(), i -> A[i][c]);
    }
    
    public float val(int r, int c) {
        return A[r][c];
    }
    
    public int width() {
        return A[0].length;
    }
    
    public int height() {
        return A.length;
    }
    
    @Override
    public String toString() {
        return Arrays.deepToString(A);
    }
    
    public static Mat identity(int size) {
        return new Mat(size, size,
            (r, c) -> r==c ? 1.0F : 0.0F);
    }
    
    public static Mat zero(int size) {
        return new Mat(size, size,
            (r, c) -> 0.0F);
    }
    
    public static Mat rows(Vec... rows) {
        return new Mat(rows.length, rows[0].dimensions(),
            (r, c) -> rows[r].A[c]);
    }
    
    public static Mat cols(Vec... cols) {
        return new Mat(cols[0].dimensions(), cols.length,
            (r, c) -> cols[c].A[r]);
    }
    
    public static Mat translate(Vec t) {
        return new Mat(t.dimensions()+1, t.dimensions()+1,
            (r, c) -> r!=c ? (c==t.dimensions() ? t.A[r] : 0.0F) : 1.0F);
     }
    
    public static Mat translate(float... t) {
        return translate(new Vec(t));
    }
    
    public static Mat rotate2(float angle) {
        
        float c = (float) Math.cos(angle);
        float s = (float) Math.sin(angle);
        
        return new Mat(new float[][] {
            {c, -s, 0.0F},
            {s, c, 0.0F},
            {0.0F, 0.0F, 1.0F}
        });
    }
    
    public static Mat rotate2(float angle, Vec point) {
        return translate(point.negate())
            .mul(rotate2(angle))
            .mul(translate(point));
    }
    
    public static Mat rotate3(float angle, Vec axis) {
        
        if(axis.dimensions() != 3)
            throw new IllegalArgumentException();
        
        Vec a = axis.normalize();
        float x = a.val(0), y = a.val(1), z = a.val(2);
        
        float c = (float) Math.cos(angle);
        float s = (float) Math.sin(angle);
        float C = 1 - c;
        
        return new Mat(new float[][] {
            {x*x*C+c, x*y*C-z*s, x*z*C+y*s, 0.0F},
            {y*x*C+z*s, y*y*C+c, y*z*C-x*s, 0.0F},
            {z*x*C-y*s, x*y*C+x*s, z*z*C+c, 0.0F},
            {0.0F, 0.0F, 0.0F, 1.0F}
        });
    }
    
    public static Mat rotate3(float angle, int axis) {
        return rotate3(angle, Vec.basis(axis, 3));
    }
    
    public static Mat rotate3(float angle, Vec axis, Vec point) {
        return translate(point.negate())
            .mul(rotate3(angle, axis))
            .mul(translate(point));
    }
    
    public static Mat rotate3(float yaw, float pitch, float roll) {
        return rotate3(yaw, Vec.basis(0, 3))
            .mul(rotate3(-pitch, Vec.basis(1, 3)))
            .mul(rotate3(roll, Vec.basis(2, 3)));
    }
    
    public static Mat scale(Vec s) {
        return new Mat(s.dimensions()+1, s.dimensions()+1,
            (r, c) -> r==c ? (r<s.dimensions() ? s.A[r] : 1.0F) : 0.0F);
    }
    
    public static Mat scale(float... s) {
        return scale(new Vec(s));
    }
    
    public static Mat projection(float fov, float aspect, float near, float far) {
        
        float ys = (float) ((1.0F / Math.tan(fov / 2.0F)) * aspect);
        float xs = ys / aspect;
        float fl = far - near;
        
        return new Mat(new float[][] {
            {xs, 0.0F, 0.0F, 0.0F},
            {0.0F, ys, 0.0F, 0.0F},
            {0.0F, 0.0F, -(far+near)/fl, -2*near*far/fl},
            {0.0F, 0.0F, -1.0F, 0.0F}
        });
    }
}