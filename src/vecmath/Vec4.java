package vecmath;

public class Vec4 {
    
    public static final Vec4 ZERO = new Vec4(0.0F, 0.0F, 0.0F, 1.0F);
    
    public static final Vec4 ONE = new Vec4(1.0F, 1.0F, 1.0F, 1.0F);
    
    public final float X, Y, Z, W;
    
    public Vec4(float x, float y, float z, float w) {
        X = x; Y = y; Z = z; W = w;
    }
    
    public Vec4(float x, float y, float z) {
        X = x; Y = y; Z = z; W = 1.0F;
    }
    
    public Vec4(Vec4 v) {
        X = v.X; Y = v.Y; Z = v.Z; W = v.W;
    }
    
    public Vec4() {
        X = 0.0F; Y = 0.0F; Z = 0.0F; W = 1.0F;
    }
    
    public Vec4 toCartesian() {
        return new Vec4(X/W, Y/W, Z/W, 1.0F);
    }
    
    public Vec4 add(Vec4 v) {
        Vec4 v1 = toCartesian(), v2 = v.toCartesian();
        return new Vec4(v1.X+v2.X, v1.Y+v2.Y, v1.Z+v2.Z, 1.0F);
    }
    
    public Vec4 sub(Vec4 v) {
        return add(v.negate());
    }
    
    public Vec4 mul(float s) {
        return new Vec4(X*s, Y*s, Z*s, W);
    }
    
    public Vec4 div(float s) {
        return mul(1.0F / s);
    }
    
    public Vec4 negate() {
        return mul(-1.0F);
    }
    
    public Vec4 normalize() {
        return mul(1.0F / length());
    }
    
    public Vec4 setLength(float length) {
        return normalize().mul(length);
    }
    
    public float length() {
        return (float) Math.sqrt((X*X + Y*Y + Z*Z) / W*W);
    }
    
    public float lengthSquared() {
        return (X*X + Y*Y + Z*Z) / W*W;
    }
    
    public float distanceTo(Vec4 v) {
        return sub(v).length();
    }
    
    public float dot(Vec4 v) {
        Vec4 v1 = toCartesian(), v2 = v.toCartesian();
        return v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z;
    }
    
    public Vec4 cross(Vec4 v) {
        
        Vec4 v1 = toCartesian(), v2 = v.toCartesian();
        
        return new Vec4(
            v1.Y*v2.Z - v1.Z*v2.Y,
            v1.Z*v2.X - v1.X*v2.Z,
            v1.X*v2.Y - v1.Y*v2.Z, 1.0F);
    }
    
    public float angle(Vec4 v) {
        return (float) Math.acos(normalize().dot(v.normalize()));
    }
}