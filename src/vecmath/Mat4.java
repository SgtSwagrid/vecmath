package vecmath;

public class Mat4 {
    
    public static final Mat4 IDENTITY = new Mat4(new float[][] {
        {1.0F, 0.0F, 0.0F, 0.0F},
        {0.0F, 1.0F, 0.0F, 0.0F},
        {0.0F, 0.0F, 1.0F, 0.0F},
        {0.0F, 0.0F, 0.0F, 1.0F}});
        
    public static final Mat4 ZERO = new Mat4(new float[][] {
        {0.0F, 0.0F, 0.0F, 0.0F},
        {0.0F, 0.0F, 0.0F, 0.0F},
        {0.0F, 0.0F, 0.0F, 0.0F},
        {0.0F, 0.0F, 0.0F, 0.0F}});
    
    private final float[][] A;
    
    public Mat4(float[][] A) { this.A = A; }
    
    public Mat4(Mat4 m) {
        
        A = new float[4][4];
        
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                A[i][j] = m.A[i][j];
            }
        }
    }
    
    public Mat4 transpose() {
        
        float[][] a = new float[4][4];
        
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                a[i][j] = A[j][i];
            }
        }
        return new Mat4(a);
    }
    
    public Mat4 add(Mat4 m) {
        
        float[][] a = new float[4][4];
        
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                a[i][j] = A[i][j] + m.A[i][j];
            }
        }
        return new Mat4(a);
    }
    
    public Mat4 sub(Mat4 m) {
        return add(m.negate());
    }
    
    public Mat4 mul(Mat4 m) {
        
        float[][] a = new float[4][4];
        
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                for(int k = 0; k < 4; k++) {
                    a[i][j] += A[i][k] * m.A[k][j];
                }
            }
        }
        return new Mat4(a);
    }
    
    public Mat4 mul(float s) {
        
        float[][] a = new float[4][4];
        
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                a[i][j] = A[i][j] * s;
            }
        }
        return new Mat4(a);
    }
    
    public Vec4 mul(Vec4 v) {
        
        float[] vec = new float[] {v.X, v.Y, v.Z, v.W};
        float[] a = new float[4];
        
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                a[i] += A[i][j] * vec[i];
            }
        }
        return new Vec4(a[0], a[1], a[2], a[3]);
    }
    
    public Mat4 div(float s) {
        return mul(1.0F / s);
    }
    
    public Mat4 negate() {
        return mul(-1.0F);
    }
    
    public float det() {
        //TODO
        return 1.0F;
    }
    
    @Override
    public String toString() {
        return A.toString();
    }
    
    public static Mat4 translate(float tx, float ty, float tz) {
        return new Mat4(new float[][] {
            {1.0F, 0.0F, 0.0F, tx},
            {0.0F, 1.0F, 0.0F, ty},
            {0.0F, 0.0F, 1.0F, tz},
            {0.0F, 0.0F, 0.0F, 1.0F}
        });
    }
    
    public static Mat4 rotate(Vec4 axis, float angle) {
        
        return null;
    }
    
    public static Mat4 scale(float sx, float sy, float sz) {
        return new Mat4(new float[][] {
            {sx, 0.0F, 0.0F, 0.0F},
            {0.0F, sy, 0.0F, 0.0F},
            {0.0F, 0.0F, sz, 0.0F},
            {0.0F, 0.0F, 0.0F, 1.0F}
        });
    }
}