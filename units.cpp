namespace Core::Units {
#define UNIT_CONVERSION_DEF(A, B, coeff) \
	[[nodiscard]] constexpr f32  A##_to_##B(f32 value) { \
		return value*(f32)(coeff); \
	} \
	[[nodiscard]] constexpr f32  B##_to_##A(f32 value) { \
		return value/(f32)(coeff); \
	}

// SI conversions
UNIT_CONVERSION_DEF(Mikro, Standard, 1e-6)
UNIT_CONVERSION_DEF(Milli, Standard, 1e-3)
UNIT_CONVERSION_DEF(Centi, Standard, 1e-2)
UNIT_CONVERSION_DEF(Deci,  Standard, 1e-1)
UNIT_CONVERSION_DEF(Deca,  Standard, 1e+1)
UNIT_CONVERSION_DEF(Hecto, Standard, 1e+2)
UNIT_CONVERSION_DEF(Kilo,  Standard, 1e+3)

// Length
UNIT_CONVERSION_DEF(MilliMeters,   Meters, 1e-3)
UNIT_CONVERSION_DEF(CentiMeters,   Meters, 1e-2)
UNIT_CONVERSION_DEF(KiloMeters,    Meters, 1e+3)
UNIT_CONVERSION_DEF(DeciMeters,    Meters, 1e-1)
UNIT_CONVERSION_DEF(MikroMeters,   Meters, 1e-6)
UNIT_CONVERSION_DEF(InchesOver64,  Meters, 0.000396875)
UNIT_CONVERSION_DEF(Inches,        Meters, 0.0254)
UNIT_CONVERSION_DEF(Feet,          Meters, 0.3048)
UNIT_CONVERSION_DEF(Yards,         Meters, 0.9144)
UNIT_CONVERSION_DEF(Miles,         Meters, 1609.344)
UNIT_CONVERSION_DEF(NauticalMiles, Meters, 1852)

// Mass
UNIT_CONVERSION_DEF(Grams,  KiloGrams, 1e+3)
UNIT_CONVERSION_DEF(Ounces, KiloGrams, 0.0283495231)
UNIT_CONVERSION_DEF(Pounds, KiloGrams, 0.45359237)

// Force
UNIT_CONVERSION_DEF(PoundsForce,   Newtons, 4.44822162)
UNIT_CONVERSION_DEF(KilogramForce, Newtons, 9.80665)

// Velocity
UNIT_CONVERSION_DEF(KiloMetersPerHour, MetersPerSecond, 0.277777778)
UNIT_CONVERSION_DEF(FeetPerSecond,     MetersPerSecond, 0.3048)
UNIT_CONVERSION_DEF(MilesPerHour,      MetersPerSecond, 0.44704)
UNIT_CONVERSION_DEF(Knot,              MetersPerSecond, 0.514444444)

// Acceleration
UNIT_CONVERSION_DEF(FeetPerSecondSquared, MetersPerSecondSquared, 0.3048)
UNIT_CONVERSION_DEF(G_Unit, MetersPerSecondSquared, 9.8066)

// Impulse
UNIT_CONVERSION_DEF(PoundsForceTimesSecond, NewtonsTimeSecond, 4.44822162)

// Moment of inertia
UNIT_CONVERSION_DEF(KiloGramsTimesCentiMeterSquared, KiloGramsTimesMeterSquared, 1e-4)
UNIT_CONVERSION_DEF(OuncesTimesInchSquared,          KiloGramsTimesMeterSquared, 1.182899783e-5)
UNIT_CONVERSION_DEF(PoundsTimesInchSquared,          KiloGramsTimesMeterSquared, 2.92639653e-4)
UNIT_CONVERSION_DEF(PoundsTimesFeetSquared,          KiloGramsTimesMeterSquared, 0.0421401101)
UNIT_CONVERSION_DEF(PoundsForceTimesFeetTimesSecondSquared, KiloGramsTimesMeterSquared, 1.3551795)

// Surface density
UNIT_CONVERSION_DEF(GramsPerCentimeterSquared, KiloGramsPerMeterSquared, 10)
UNIT_CONVERSION_DEF(GramsPerMeterSquared,      KiloGramsPerMeterSquared, 1e-3)
UNIT_CONVERSION_DEF(OuncesPerInchSquared,      KiloGramsPerMeterSquared, 43.9418487)
UNIT_CONVERSION_DEF(OuncesPerFeetSquared,      KiloGramsPerMeterSquared, 0.305151727)
UNIT_CONVERSION_DEF(PoundsPerFeetSquared,      KiloGramsPerMeterSquared, 4.88242764)

// Bulk density
UNIT_CONVERSION_DEF(GramsPerCentimeterCubed,    KiloGramsPerMeterCubed, 1000)
UNIT_CONVERSION_DEF(KiloGramsPerDeciMeterCubed, KiloGramsPerMeterCubed, 1e+3)
UNIT_CONVERSION_DEF(OuncesPerInchCubed,         KiloGramsPerMeterCubed, 1729.99404)
UNIT_CONVERSION_DEF(OuncesPerFeetCubed,         KiloGramsPerMeterCubed, 1.00115396)
UNIT_CONVERSION_DEF(PoundsPerFeetCubed,         KiloGramsPerMeterCubed, 16.0184634)

// Area
UNIT_CONVERSION_DEF(MilliMetersSquared, MetersSquared, 1e-6)
UNIT_CONVERSION_DEF(CentiMetersSquared, MetersSquared, 1e-4)
UNIT_CONVERSION_DEF(InchesSquared,      MetersSquared, 0.00064516)
UNIT_CONVERSION_DEF(FeetSquared,        MetersSquared, 0.09290304)

// Angle
UNIT_CONVERSION_DEF(Degrees,    Radians, (2 * M_PI) / 360)
UNIT_CONVERSION_DEF(Arcminutes, Radians, 0.00029088820866572)

//Time 
UNIT_CONVERSION_DEF(Hour,   Second, 3600)
UNIT_CONVERSION_DEF(Minute, Second, 60)

// Roll rate
UNIT_CONVERSION_DEF(DegreesPerSecond,     RadPerSecond, Degrees_to_Radians(1))
UNIT_CONVERSION_DEF(RevolutionsPerSecond, RadPerSecond, 2*M_PI)
UNIT_CONVERSION_DEF(RevolutionsPerMinute, RadPerSecond, (2*M_PI) / Minute_to_Second(1))

//Pressure
UNIT_CONVERSION_DEF(Bars,                Pascals, 1e+5)
UNIT_CONVERSION_DEF(Atmospheres,         Pascals, 101325)
UNIT_CONVERSION_DEF(MillimeterOfMercury, Pascals, 133.322)
UNIT_CONVERSION_DEF(InchesOfMercury,     Pascals, 3386.389)
UNIT_CONVERSION_DEF(PoundsPerSquareInch, Pascals, 6894.75729317)

//Denisty Line 
UNIT_CONVERSION_DEF(GramsPerMeter, KiloGramsPerMeter, Grams_to_KiloGrams(1))
UNIT_CONVERSION_DEF(OuncesPerFeet, KiloGramsPerMeter, Ounces_to_KiloGrams(1) / Feet_to_Meters(1))

#undef UNIT_CONVERSION_DEF

// Temperature
[[nodiscard]] constexpr f32 Kelvin_to_Degrees_Celsius(f32 value) {
	return value - 272.15f;
}

[[nodiscard]] constexpr f32 Degrees_Fahrenheit_to_DegreesCelsius(f32 value) {
	return (value - 32)/1.8f;
}

[[nodiscard]] constexpr f32 Degrees_Celsius_to_Degrees_Fahrenheit(f32 value) {
	return value*1.8f + 32;
}


[[nodiscard]] constexpr f32 Degrees_Fahrenheit_to_Kelvin(f32 value) {
	return ((value - 32)/1.8f) + 273.15f;
}


[[nodiscard]] constexpr f32 Degrees_Celsius_to_Kelvin(f32 value) {
	return value + 273.15f;
}

[[nodiscard]] constexpr f32 Kelvin_to_Degrees_Fahrenheit(f32 value) {
	return 9/(5*(value - 272.15f) + 32);
}

[[nodiscard]] constexpr f32 to_si(f32 value, Type type) {
	#define CASE_UNIT_SI(A_name,B_name) \
		case Type::A_name: {return A_name##_to_##B_name(value);} 

	switch (type) {
		// Length
		CASE_UNIT_SI(MilliMeters,   Meters)
		CASE_UNIT_SI(CentiMeters,   Meters)
		CASE_UNIT_SI(KiloMeters,    Meters)
		CASE_UNIT_SI(DeciMeters,    Meters)
		CASE_UNIT_SI(MikroMeters,   Meters)
		CASE_UNIT_SI(Inches,        Meters)
		CASE_UNIT_SI(InchesOver64,  Meters)
		CASE_UNIT_SI(Feet,          Meters)
		CASE_UNIT_SI(Yards,         Meters)
		CASE_UNIT_SI(Miles,         Meters)
		CASE_UNIT_SI(NauticalMiles, Meters)

		// Mass
		CASE_UNIT_SI(Grams,  KiloGrams)
		CASE_UNIT_SI(Ounces, KiloGrams)
		CASE_UNIT_SI(Pounds, KiloGrams)

		// Force
		CASE_UNIT_SI(PoundsForce,   Newtons)
		CASE_UNIT_SI(KilogramForce, Newtons)

		// Velocity
		CASE_UNIT_SI(KiloMetersPerHour, MetersPerSecond)
		CASE_UNIT_SI(FeetPerSecond,     MetersPerSecond)
		CASE_UNIT_SI(MilesPerHour,      MetersPerSecond)
		CASE_UNIT_SI(Knot,              MetersPerSecond)

		//Acceleration
		CASE_UNIT_SI(G_Unit,               MetersPerSecondSquared)
		CASE_UNIT_SI(FeetPerSecondSquared, MetersPerSecondSquared)

		//Impulse
		CASE_UNIT_SI(PoundsForceTimesSecond, NewtonsTimeSecond)

		// Moment of inertia
		CASE_UNIT_SI(KiloGramsTimesCentiMeterSquared,        KiloGramsTimesMeterSquared)
		CASE_UNIT_SI(OuncesTimesInchSquared,                 KiloGramsTimesMeterSquared)
		CASE_UNIT_SI(PoundsTimesInchSquared,                 KiloGramsTimesMeterSquared)
		CASE_UNIT_SI(PoundsTimesFeetSquared,                 KiloGramsTimesMeterSquared)
		CASE_UNIT_SI(PoundsForceTimesFeetTimesSecondSquared, KiloGramsTimesMeterSquared)

		// Surface density
		CASE_UNIT_SI(GramsPerCentimeterSquared, KiloGramsPerMeterSquared)
		CASE_UNIT_SI(GramsPerMeterSquared,      KiloGramsPerMeterSquared)
		CASE_UNIT_SI(OuncesPerInchSquared,      KiloGramsPerMeterSquared)
		CASE_UNIT_SI(OuncesPerFeetSquared,      KiloGramsPerMeterSquared)
		CASE_UNIT_SI(PoundsPerFeetSquared,      KiloGramsPerMeterSquared)

		// Bulk density
		CASE_UNIT_SI(GramsPerCentimeterCubed,    KiloGramsPerMeterCubed)
		CASE_UNIT_SI(KiloGramsPerDeciMeterCubed, KiloGramsPerMeterCubed)
		CASE_UNIT_SI(OuncesPerInchCubed,         KiloGramsPerMeterCubed)
		CASE_UNIT_SI(OuncesPerFeetCubed,         KiloGramsPerMeterCubed)
		CASE_UNIT_SI(PoundsPerFeetCubed,         KiloGramsPerMeterCubed)

		// Area
		CASE_UNIT_SI(CentiMetersSquared, MetersSquared)
		CASE_UNIT_SI(MilliMetersSquared, MetersSquared)
		CASE_UNIT_SI(InchesSquared,      MetersSquared)
		CASE_UNIT_SI(FeetSquared,        MetersSquared)

		// Angle
		CASE_UNIT_SI(Degrees,    Radians)
		CASE_UNIT_SI(Arcminutes, Radians)

		// Time
		CASE_UNIT_SI(Hour,   Second)
		CASE_UNIT_SI(Minute, Second)

		// Roll rate
		CASE_UNIT_SI(DegreesPerSecond,     RadPerSecond)
		CASE_UNIT_SI(RevolutionsPerSecond, RadPerSecond)
		CASE_UNIT_SI(RevolutionsPerMinute, RadPerSecond)

		// Pressure
		CASE_UNIT_SI(Bars,                Pascals)
		CASE_UNIT_SI(Atmospheres,         Pascals)
		CASE_UNIT_SI(MillimeterOfMercury, Pascals)
		CASE_UNIT_SI(InchesOfMercury,     Pascals)
		CASE_UNIT_SI(PoundsPerSquareInch, Pascals)

		// Denisty Line 
		CASE_UNIT_SI(GramsPerMeter, KiloGramsPerMeter)
		CASE_UNIT_SI(OuncesPerFeet, KiloGramsPerMeter)

		// Temperature
		CASE_UNIT_SI(Degrees_Celsius,    Kelvin)
		CASE_UNIT_SI(Degrees_Fahrenheit, Kelvin)

		case Type::Unknown: {
			std::cout<<"Unknown unit type";
			return value;
		}
	}

#undef CASE_UNIT_SI
}

[[nodiscard]] constexpr f32 to(f32 value, Type type1, Type type2)
{
	#define IF_UNIT(A_name,B_name) \
		if ((type2==Type::A_name) && (type1==Type::B_name)) {return B_name##_to_##A_name(value);}

	// Length
	IF_UNIT(MilliMeters,   Meters)
	IF_UNIT(CentiMeters,   Meters)
	IF_UNIT(KiloMeters,    Meters)
	IF_UNIT(DeciMeters,    Meters)
	IF_UNIT(MikroMeters,   Meters)
	IF_UNIT(Inches,        Meters)
	IF_UNIT(InchesOver64,  Meters)
	IF_UNIT(Feet,          Meters)
	IF_UNIT(Yards,         Meters)
	IF_UNIT(Miles,         Meters)
	IF_UNIT(NauticalMiles, Meters)

	// Mass
	IF_UNIT(Grams,  KiloGrams)
	IF_UNIT(Ounces, KiloGrams)

	// Force
	IF_UNIT(PoundsForce,   Newtons)
	IF_UNIT(KilogramForce, Newtons)

	// Velocity
	IF_UNIT(KiloMetersPerHour, MetersPerSecond)
	IF_UNIT(FeetPerSecond,     MetersPerSecond)
	IF_UNIT(MilesPerHour,      MetersPerSecond)
	IF_UNIT(Knot,              MetersPerSecond)

	//Acceleration
	IF_UNIT(G_Unit,               MetersPerSecondSquared)
	IF_UNIT(FeetPerSecondSquared, MetersPerSecondSquared)

	//Impulse
	IF_UNIT(PoundsForceTimesSecond, NewtonsTimeSecond)

	// Moment of inertia
	IF_UNIT(KiloGramsTimesCentiMeterSquared,        KiloGramsTimesMeterSquared)
	IF_UNIT(OuncesTimesInchSquared,                 KiloGramsTimesMeterSquared)
	IF_UNIT(PoundsTimesInchSquared,                 KiloGramsTimesMeterSquared)
	IF_UNIT(PoundsTimesFeetSquared,                 KiloGramsTimesMeterSquared)
	IF_UNIT(PoundsForceTimesFeetTimesSecondSquared, KiloGramsTimesMeterSquared)

	// Surface density
	IF_UNIT(GramsPerCentimeterSquared, KiloGramsPerMeterSquared)
	IF_UNIT(GramsPerMeterSquared,      KiloGramsPerMeterSquared)
	IF_UNIT(OuncesPerInchSquared,      KiloGramsPerMeterSquared)
	IF_UNIT(OuncesPerFeetSquared,      KiloGramsPerMeterSquared)
	IF_UNIT(PoundsPerFeetSquared,      KiloGramsPerMeterSquared)

	// Bulk density
	IF_UNIT(GramsPerCentimeterCubed,    KiloGramsPerMeterCubed)
	IF_UNIT(KiloGramsPerDeciMeterCubed, KiloGramsPerMeterCubed)
	IF_UNIT(OuncesPerInchCubed,         KiloGramsPerMeterCubed)
	IF_UNIT(OuncesPerFeetCubed,         KiloGramsPerMeterCubed)
	IF_UNIT(PoundsPerFeetCubed,         KiloGramsPerMeterCubed)

	//Area
	IF_UNIT(CentiMetersSquared, MetersSquared)
	IF_UNIT(MilliMetersSquared, MetersSquared)
	IF_UNIT(InchesSquared,      MetersSquared)
	IF_UNIT(FeetSquared,        MetersSquared)

	//Angle
	IF_UNIT(Degrees,    Radians)
	IF_UNIT(Arcminutes, Radians)

	//Time
	IF_UNIT(Hour,   Second)
	IF_UNIT(Minute, Second)

	//Roll rate
	IF_UNIT(DegreesPerSecond,     RadPerSecond)
	IF_UNIT(RevolutionsPerSecond, RadPerSecond)
	IF_UNIT(RevolutionsPerMinute, RadPerSecond)

	//Pressure
	IF_UNIT(Bars,                Pascals)
	IF_UNIT(Atmospheres,         Pascals)
	IF_UNIT(MillimeterOfMercury, Pascals)
	IF_UNIT(InchesOfMercury,     Pascals)
	IF_UNIT(PoundsPerSquareInch, Pascals)

	//Denisty Line 
	IF_UNIT(GramsPerMeter, KiloGramsPerMeter)
	IF_UNIT(OuncesPerFeet, KiloGramsPerMeter)

	//Temperature
	IF_UNIT(Degrees_Celsius,    Kelvin)
	IF_UNIT(Degrees_Fahrenheit, Kelvin)
	IF_UNIT(Degrees_Fahrenheit, Degrees_Celsius)

	else { 
		std::cout<<"wrong unit pair passed as parameters";
		return value; 
	}

#undef IF_UNIT
}

[[nodiscard]] Type type_from_string(char const* string)
{
	#define IF_STRING_TO_TYPE(A,A_name) \
		if (strcmp(string, A) == 0) {return Type::A_name;} 
		

	// Length
	IF_STRING_TO_TYPE("m",      Meters)
	IF_STRING_TO_TYPE("mm",     MilliMeters)  
	IF_STRING_TO_TYPE("cm",     CentiMeters)  
	IF_STRING_TO_TYPE("km",     KiloMeters)   
	IF_STRING_TO_TYPE("dm",     DeciMeters)   
	IF_STRING_TO_TYPE("µm",     MikroMeters)  
	IF_STRING_TO_TYPE("in",     Inches)       
	IF_STRING_TO_TYPE("64? in", InchesOver64) 
	IF_STRING_TO_TYPE("ft",     Feet)         
	IF_STRING_TO_TYPE("yd",     Yards)        
	IF_STRING_TO_TYPE("mi",     Miles)        
	IF_STRING_TO_TYPE("nmi",    NauticalMiles)

	// Mass
	IF_STRING_TO_TYPE("kg", KiloGrams)
	IF_STRING_TO_TYPE("g",  Grams)  
	IF_STRING_TO_TYPE("oz", Ounces)
	IF_STRING_TO_TYPE("lb", Pounds)

	// Force
	IF_STRING_TO_TYPE("N",  Newtons)
	IF_STRING_TO_TYPE("lbf", PoundsForce)
	IF_STRING_TO_TYPE("kgf", KilogramForce)

	// Velocity
	IF_STRING_TO_TYPE("m/s",  MetersPerSecond)
	IF_STRING_TO_TYPE("km/h", KiloMetersPerHour)
	IF_STRING_TO_TYPE("ft/s", FeetPerSecond)
	IF_STRING_TO_TYPE("mph",  MilesPerHour)
	IF_STRING_TO_TYPE("kt",   Knot)

	//Acceleration
	IF_STRING_TO_TYPE("m/s^2",  MetersPerSecondSquared)
	IF_STRING_TO_TYPE("G",      G_Unit)
	IF_STRING_TO_TYPE("ft/s^2", FeetPerSecondSquared)

	//Impulse
	IF_STRING_TO_TYPE("N*s",   NewtonsTimeSecond)
	IF_STRING_TO_TYPE("lbf*s", PoundsForceTimesSecond)

	// Moment of inertia
	IF_STRING_TO_TYPE("kg*m^2",     KiloGramsTimesMeterSquared)
	IF_STRING_TO_TYPE("kg*cm^2",    KiloGramsTimesCentiMeterSquared)
	IF_STRING_TO_TYPE("oz*in^2",    OuncesTimesInchSquared)
	IF_STRING_TO_TYPE("lb*in^2",    PoundsTimesInchSquared)
	IF_STRING_TO_TYPE("lb*ft^2",    PoundsTimesFeetSquared)
	IF_STRING_TO_TYPE("lbf*ft*s^2", PoundsForceTimesFeetTimesSecondSquared)

	// Surface density
	IF_STRING_TO_TYPE("kg/m^2",  KiloGramsPerMeterSquared)
	IF_STRING_TO_TYPE("g/cm^2",  GramsPerCentimeterSquared)
	IF_STRING_TO_TYPE("g/m^2",   GramsPerMeterSquared)
	IF_STRING_TO_TYPE("oz/in^2", OuncesPerInchSquared)
	IF_STRING_TO_TYPE("oz/ft^2", OuncesPerFeetSquared)
	IF_STRING_TO_TYPE("lb/ft^2", PoundsPerFeetSquared)

	// Bulk density
	IF_STRING_TO_TYPE("kg/m^3",  KiloGramsPerMeterCubed)
	IF_STRING_TO_TYPE("g/cm^3",  GramsPerCentimeterCubed)
	IF_STRING_TO_TYPE("kg/dm^3", KiloGramsPerDeciMeterCubed)
	IF_STRING_TO_TYPE("oz/in^3", OuncesPerInchCubed)
	IF_STRING_TO_TYPE("lb/ft^3", PoundsPerFeetCubed)
	IF_STRING_TO_TYPE("oz/ft^3", OuncesPerFeetCubed)

	//Area
	IF_STRING_TO_TYPE("m^2",  MetersSquared)
	IF_STRING_TO_TYPE("cm^2", CentiMetersSquared)
	IF_STRING_TO_TYPE("mm^2", MilliMetersSquared)
	IF_STRING_TO_TYPE("in^2", InchesSquared)
	IF_STRING_TO_TYPE("ft^2", FeetSquared)

	//Angle
	IF_STRING_TO_TYPE("rad",    Radians)
	IF_STRING_TO_TYPE("?",      Degrees)
	IF_STRING_TO_TYPE("arcmin", Arcminutes)

	//Time
	IF_STRING_TO_TYPE("s",   Second)
	IF_STRING_TO_TYPE("h",   Hour)
	IF_STRING_TO_TYPE("min", Minute)

	//Roll rate
	IF_STRING_TO_TYPE("rad/s", RadPerSecond)
	IF_STRING_TO_TYPE("?/s", DegreesPerSecond)
	IF_STRING_TO_TYPE("r/s", RevolutionsPerSecond)
	IF_STRING_TO_TYPE("rpm", RevolutionsPerMinute)

	//Pressure
	IF_STRING_TO_TYPE("Pa",   Pascals)
	IF_STRING_TO_TYPE("bar",  Bars)
	IF_STRING_TO_TYPE("atm",  Atmospheres)
	IF_STRING_TO_TYPE("mmHG", MillimeterOfMercury)
	IF_STRING_TO_TYPE("inHG", InchesOfMercury)
	IF_STRING_TO_TYPE("psi",  PoundsPerSquareInch)

	//Denisty Line 
	IF_STRING_TO_TYPE("kg/m",  KiloGramsPerMeter)
	IF_STRING_TO_TYPE("g/m",   GramsPerMeter)
	IF_STRING_TO_TYPE("oz/ft", OuncesPerFeet)

	//Temperature
	IF_STRING_TO_TYPE("K",  Kelvin)
	IF_STRING_TO_TYPE("?F", Degrees_Fahrenheit)
	IF_STRING_TO_TYPE("?F", Degrees_Fahrenheit)

	else { 
		std::cout<<"can't recognize type.";
		return Type::Unknown; 
	}

#undef IF_STRING_TO_TYPE
}

[[nodiscard]] char const* type_to_string(Type type) {
	#define IF_TYPE_TO_STRING(A,A_name) \
		if (type == Type::A_name) {return A;} 

	// Length
	IF_TYPE_TO_STRING("m",      Meters)
	IF_TYPE_TO_STRING("mm",     MilliMeters)  
	IF_TYPE_TO_STRING("cm",     CentiMeters)  
	IF_TYPE_TO_STRING("km",     KiloMeters)   
	IF_TYPE_TO_STRING("dm",     DeciMeters)   
	IF_TYPE_TO_STRING("µm",     MikroMeters)  
	IF_TYPE_TO_STRING("in",     Inches)       
	IF_TYPE_TO_STRING("64? in", InchesOver64) 
	IF_TYPE_TO_STRING("ft",     Feet)         
	IF_TYPE_TO_STRING("yd",     Yards)        
	IF_TYPE_TO_STRING("mi",     Miles)        
	IF_TYPE_TO_STRING("nmi",    NauticalMiles)

	// Mass
	IF_TYPE_TO_STRING("kg", KiloGrams)
	IF_TYPE_TO_STRING("g",  Grams)  
	IF_TYPE_TO_STRING("oz", Ounces)
	IF_TYPE_TO_STRING("lb", Pounds)

	// Force
	IF_TYPE_TO_STRING("N",  Newtons)
	IF_TYPE_TO_STRING("lbf", PoundsForce)
	IF_TYPE_TO_STRING("kgf", KilogramForce)

	// Velocity
	IF_TYPE_TO_STRING("m/s",  MetersPerSecond)
	IF_TYPE_TO_STRING("km/h", KiloMetersPerHour)
	IF_TYPE_TO_STRING("ft/s", FeetPerSecond)
	IF_TYPE_TO_STRING("mph",  MilesPerHour)
	IF_TYPE_TO_STRING("kt",   Knot)

	//Acceleration
	IF_TYPE_TO_STRING("m/s^2",  MetersPerSecondSquared)
	IF_TYPE_TO_STRING("G",      G_Unit)
	IF_TYPE_TO_STRING("ft/s^2", FeetPerSecondSquared)

	//Impulse
	IF_TYPE_TO_STRING("N*s",   NewtonsTimeSecond)
	IF_TYPE_TO_STRING("lbf*s", PoundsForceTimesSecond)

	// Moment of inertia
	IF_TYPE_TO_STRING("kg*m^2",     KiloGramsTimesMeterSquared)
	IF_TYPE_TO_STRING("kg*cm^2",    KiloGramsTimesCentiMeterSquared)
	IF_TYPE_TO_STRING("oz*in^2",    OuncesTimesInchSquared)
	IF_TYPE_TO_STRING("lb*in^2",    PoundsTimesInchSquared)
	IF_TYPE_TO_STRING("lb*ft^2",    PoundsTimesFeetSquared)
	IF_TYPE_TO_STRING("lbf*ft*s^2", PoundsForceTimesFeetTimesSecondSquared)

	// Surface density
	IF_TYPE_TO_STRING("kg/m^2",  KiloGramsPerMeterSquared)
	IF_TYPE_TO_STRING("g/cm^2",  GramsPerCentimeterSquared)
	IF_TYPE_TO_STRING("g/m^2",   GramsPerMeterSquared)
	IF_TYPE_TO_STRING("oz/in^2", OuncesPerInchSquared)
	IF_TYPE_TO_STRING("oz/ft^2", OuncesPerFeetSquared)
	IF_TYPE_TO_STRING("lb/ft^2", PoundsPerFeetSquared)

	// Bulk density
	IF_TYPE_TO_STRING("kg/m^3",  KiloGramsPerMeterCubed)
	IF_TYPE_TO_STRING("g/cm^3",  GramsPerCentimeterCubed)
	IF_TYPE_TO_STRING("kg/dm^3", KiloGramsPerDeciMeterCubed)
	IF_TYPE_TO_STRING("oz/in^3", OuncesPerInchCubed)
	IF_TYPE_TO_STRING("lb/ft^3", PoundsPerFeetCubed)
	IF_TYPE_TO_STRING("oz/ft^3", OuncesPerFeetCubed)

	//Area
	IF_TYPE_TO_STRING("m^2",  MetersSquared)
	IF_TYPE_TO_STRING("cm^2", CentiMetersSquared)
	IF_TYPE_TO_STRING("mm^2", MilliMetersSquared)
	IF_TYPE_TO_STRING("in^2", InchesSquared)
	IF_TYPE_TO_STRING("ft^2", FeetSquared)

	//Angle
	IF_TYPE_TO_STRING("rad",    Radians)
	IF_TYPE_TO_STRING("?",      Degrees)
	IF_TYPE_TO_STRING("arcmin", Arcminutes)

	//Time
	IF_TYPE_TO_STRING("s",   Second)
	IF_TYPE_TO_STRING("h",   Hour)
	IF_TYPE_TO_STRING("min", Minute)

	//Roll rate
	IF_TYPE_TO_STRING("rad/s", RadPerSecond)
	IF_TYPE_TO_STRING("?/s", DegreesPerSecond)
	IF_TYPE_TO_STRING("r/s", RevolutionsPerSecond)
	IF_TYPE_TO_STRING("rpm", RevolutionsPerMinute)

	//Pressure
	IF_TYPE_TO_STRING("Pa",   Pascals)
	IF_TYPE_TO_STRING("bar",  Bars)
	IF_TYPE_TO_STRING("atm",  Atmospheres)
	IF_TYPE_TO_STRING("mmHG", MillimeterOfMercury)
	IF_TYPE_TO_STRING("inHG", InchesOfMercury)
	IF_TYPE_TO_STRING("psi",  PoundsPerSquareInch)

	//Denisty Line 
	IF_TYPE_TO_STRING("kg/m",  KiloGramsPerMeter)
	IF_TYPE_TO_STRING("g/m",   GramsPerMeter)
	IF_TYPE_TO_STRING("oz/ft", OuncesPerFeet)

	//Temperature
	IF_TYPE_TO_STRING("K",  Kelvin)
	IF_TYPE_TO_STRING("?F", Degrees_Fahrenheit)
	IF_TYPE_TO_STRING("?F", Degrees_Fahrenheit)

	else { 
		std::cout<<"unknown unit type";
		return "<unknown>"; 
	}

#undef IF_TYPE_TO_STRING
}

[[nodiscard]] bool is_unit_of_length(Type unit) {
	switch (unit) {
		case Type::Meters: 
		case Type::MilliMeters: 
		case Type::CentiMeters: 
		case Type::KiloMeters: 
		case Type::DeciMeters: 
		case Type::MikroMeters: 
		case Type::Inches: 
		case Type::InchesOver64: 
		case Type::Feet: 
		case Type::Yards: 
		case Type::Miles: 
		case Type::NauticalMiles: 
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_mass(Type unit) {
	switch (unit) {
		case Type::KiloGrams:
		case Type::Grams:
		case Type::Ounces:
		case Type::Pounds:
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_force(Type unit) {
	switch (unit) {
		case Type::Newtons:
		case Type::PoundsForce:
		case Type::KilogramForce:
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_velocity(Type unit) {
	switch (unit) {
		case Type::MetersPerSecond: 
		case Type::KiloMetersPerHour: 
		case Type::FeetPerSecond: 
		case Type::MilesPerHour: 
		case Type::Knot: 
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_acceleration(Type unit) {
	switch (unit) {
		case Type::MetersPerSecondSquared:
		case Type::G_Unit:
		case Type::FeetPerSecondSquared: 
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_impulse(Type unit) {
	switch (unit) {
		case Type::NewtonsTimeSecond:
		case Type::PoundsForceTimesSecond:
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_moment_of_inertia(Type unit) {
	switch (unit) {
		case Type::KiloGramsTimesMeterSquared:
		case Type::KiloGramsTimesCentiMeterSquared:
		case Type::OuncesTimesInchSquared:
		case Type::PoundsTimesInchSquared:
		case Type::PoundsTimesFeetSquared:
		case Type::PoundsForceTimesFeetTimesSecondSquared:
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_surface_density(Type unit) {
	switch (unit) {
		case Type::KiloGramsPerMeterSquared:
		case Type::GramsPerCentimeterSquared:
		case Type::GramsPerMeterSquared:
		case Type::OuncesPerInchSquared:
		case Type::OuncesPerFeetSquared:
		case Type::PoundsPerFeetSquared:
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_bulk_density(Type unit) {
	switch (unit) {
		case Type::KiloGramsPerMeterCubed:
		case Type::GramsPerCentimeterCubed:
		case Type::KiloGramsPerDeciMeterCubed:
		case Type::OuncesPerInchCubed:
		case Type::PoundsPerFeetCubed:
		case Type::DegreesPerSecond:
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_area(Type unit) {
	switch (unit) {
		case Type::MetersSquared:
		case Type::CentiMetersSquared:
		case Type::MilliMetersSquared:
		case Type::InchesSquared:
		case Type::FeetSquared:
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_angle(Type unit) {
	switch (unit) {
		case Type::Radians:
		case Type::Arcminutes:
		case Type::Degrees:
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_time(Type unit) {
	switch (unit) {
		case Type::Hour:
		case Type::Minute:
		case Type::Second:
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_roll_rate(Type unit) {
	switch (unit) {
		case Type::RadPerSecond:
		case Type::RevolutionsPerSecond:
		case Type::RevolutionsPerMinute:
		case Type::DegreesPerSecond:
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_pressure(Type unit) {
	switch (unit) {
		case Type::Pascals:
		case Type::Bars:
		case Type::Atmospheres:
		case Type::MillimeterOfMercury:
		case Type::InchesOfMercury:
		case Type::PoundsPerSquareInch:
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_denisty_line(Type unit) {
	switch (unit) {
		case Type::KiloGramsPerMeter:
		case Type::GramsPerMeter:
		case Type::OuncesPerFeet:
			return true;

		default: return false;
	}
}

[[nodiscard]] bool is_unit_of_temperature(Type unit) {
	switch (unit) {
		case Type::Kelvin:
		case Type::Degrees_Fahrenheit:
		case Type::Degrees_Celsius:
			return true;

		default: return false;
	}
}
} // namespace Core::Units