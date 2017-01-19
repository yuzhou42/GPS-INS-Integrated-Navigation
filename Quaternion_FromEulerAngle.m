function [w,x,y,z]=Quaternion_FromEulerAngle( ea_fRoll,ea_fPitch,ea_fYaw)         %%欧拉转四元

fCosHRoll = cos(ea_fRoll * 0.5);
fSinHRoll = sin(ea_fRoll * 0.5);
fCosHPitch = cos(ea_fPitch * 0.5);
fSinHPitch = sin(ea_fPitch * 0.5);
fCosHYaw = cos(ea_fYaw * 0.5);
fSinHYaw = sin(ea_fYaw * 0.5);

w = fCosHRoll*fCosHPitch*fCosHYaw+fSinHRoll*fSinHPitch*fSinHYaw;
x = fCosHPitch * fSinHRoll * fCosHYaw - fCosHRoll * fSinHPitch * fSinHYaw;
y = fCosHRoll * fCosHYaw  * fSinHPitch + fSinHRoll * fCosHPitch * fSinHYaw;
z = fCosHRoll * fCosHPitch * fSinHYaw - fCosHYaw  * fSinHPitch * fSinHRoll;
