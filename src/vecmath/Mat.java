package vecmath;

import java.io.Serializable;
import java.util.Arrays;
import java.util.function.BiFunction;

/**
 * Arbitrary-dimensional immutable matrix type for computer graphics.
 * Used to represent general affine and projective transformations.
 * @author Alec Dorrington
 */
public class Mat implements Serializable {
    
    private static final long serialVersionUID = 3888722482803464661L;
    
    /** The contents of the matrix. */
    final float[][] A;
    
    /**
     * Construct a new matrix.
     * @param a the contents of the matrix, by row then column.
     */
    public Mat(float[][] a) {
        this(a.length, a[0].length, (r, c) -> a[r][c]);
    }
    
    /**
     * Construct a new matrix using a generator. The generator
     * accepts a row and column and should return a corresponding value.
     * @param rows the number of rows in the matrix.
     * @param width the number of columns in the matrix.
     * @param generator the function for generating the matrix's contents.
     */
    public Mat(int rows, int cols,
            BiFunction<Integer, Integer, Float> generator) {
        
        A = new float[rows][cols];
        
        for(int r = 0; r < rows; r++) {
            for(int c = 0; c < cols; c++) {
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
    
    public float trace() {
        
        if(height() != width())
            throw new IllegalArgumentException();
        
        float trace = 0.0F;
        
        for(int i = 0; i < height(); i++) {
            trace += A[i][i];
        }
        return trace;
    }
    
    public Mat changeBasis(Mat basis) {
        return basis.invert().mul(this);
    }
    
    public Vec getTranslation() {
        return col(width()-1).toCartesian();
    }
    
    public Mat getRotationMatrix() {
        return scale(getScale()).invert()
            .mul(translate(getTranslation().negate()))
            .mul(this);
    }
    
    public float getRotation2() {
        return (float) Math.atan2(getRotationMatrix().A[1][0],
            getRotationMatrix().A[0][0]);
    }
    
    public Quat getRotation3() {
        
        Mat r = getRotationMatrix();
        
        return new Quat(new float[] {
            (float) Math.sqrt(1.0F + r.A[0][0] + r.A[1][1] + r.A[2][2])
                * 0.5F,
            (float) Math.sqrt(1.0F + r.A[0][0] - r.A[1][1] - r.A[2][2])
                * 0.5F * Math.signum(r.A[2][1] - r.A[1][2]),
            (float) Math.sqrt(1.0F - r.A[0][0] + r.A[1][1] - r.A[2][2])
                * 0.5F * Math.signum(r.A[0][2] - r.A[2][0]),
            (float) Math.sqrt(1.0F - r.A[0][0] - r.A[1][1] + r.A[2][2])
                * 0.5F * Math.signum(r.A[1][0] - r.A[0][1]),
        }).normalize();
    }
    
    public Vec getScale() {
        return new Vec(width()-1, i -> col(i).length());
    }
    
    public Mat interpolate2(Mat m, float t) {
        
        Vec t1 = getTranslation(), t2 = m.getTranslation();
        float r1 = getRotation2(), r2 = m.getRotation2();
        Vec s1 = getScale(), s2 = m.getScale();
        
        float angle = r2 - r1;
        if(angle > Math.PI) angle -= 2.0F * Math.PI;
        return transform2(t1.lerp(t2, t), angle*t, s1.lerp(s2, t));
    }
    
    public Mat interpolate3(Mat m, float t) {
        
        Vec t1 = getTranslation(), t2 = m.getTranslation();
        Quat r1 = getRotation3(), r2 = m.getRotation3();
        Vec s1 = getScale(), s2 = m.getScale();
        
        return transform3(t1.lerp(t2, t), r1.slerp(r2, t), s1.lerp(s2, t));
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
    
    public static Mat rotate2(float a) {
        
        float c = (float) Math.cos(a);
        float s = (float) Math.sin(a);
        
        return new Mat(new float[][] {
            {c, -s, 0.0F},
            {s, c, 0.0F},
            {0.0F, 0.0F, 1.0F}
        });
    }
    
    public static Mat rotate2(float a, Vec o) {
        return translate(o)
            .mul(rotate2(a))
            .mul(translate(o.negate()));
    }
    
    public static Mat rotate3(Quat q) {
        
        q = q.normalize();
        float a = q.A[0], b = q.A[1], c = q.A[2], d = q.A[3];
        
        return new Mat(new float[][] {
            {a*a + b*b - c*c - d*d, 2*b*c - 2*a*d, 2*b*d + 2*a*c, 0.0F},
            {2*b*c + 2*a*d, a*a - b*b + c*c - d*d, 2*c*d - 2*a*b, 0.0F},
            {2*b*d - 2*a*c, 2*c*d + 2*a*b, a*a - b*b - c*c + d*d, 0.0F},
            {0.0F, 0.0F, 0.0F, 1.0F}
        });
    }
    
    public static Mat rotate3(Quat q, Vec o) {
        return translate(o)
            .mul(rotate3(q))
            .mul(translate(o.negate()));
    }
    
    public static Mat scale(Vec s) {
        return new Mat(s.dimensions()+1, s.dimensions()+1,
            (r, c) -> r==c ? (r<s.dimensions() ? s.A[r] : 1.0F) : 0.0F);
    }
    
    public static Mat scale(float... s) {
        return scale(new Vec(s));
    }
    
    public static Mat transform2(Vec translation, float rotation, Vec scale) {
        return translate(translation).mul(rotate2(rotation)).mul(scale(scale));
    }
    
    public static Mat transform2(Vec translation, float rotation) {
        return translate(translation).mul(rotate2(rotation));
    }
    
    public static Mat transform3(Vec translation, Quat rotation, Vec scale) {
        return translate(translation).mul(rotate3(rotation)).mul(scale(scale));
    }
    
    public static Mat transform3(Vec translation, Quat rotation) {
        return translate(translation).mul(rotate3(rotation));
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