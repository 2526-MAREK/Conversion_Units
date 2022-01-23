#include <iostream>
namespace Core::Units {
namespace Constants {
	static constexpr auto THE_CONSTANT_OF_GRAVITY = 9.80665f; // m/s^2, earth acceleration
};

enum class Type {
	Kilo,
	Hecto,
	Deca,
	Deci,
	Centi,
	Mili,
	Mikro,
	Meters,
	MilliMeters,
	CentiMeters,
	KiloMeters,
	DeciMeters,
	MikroMeters,
	Grams,
	InchesOver64,
	Inches,
	Feet,
	Yards,
	Miles,
	NauticalMiles,
	KiloGrams,
	Ounces, //oz
	Pounds, //lb
	PoundsForce,
	Newtons,
	KilogramForce,
	MetersPerSecond,
	KiloMetersPerHour,
	FeetPerSecond,
	MilesPerHour,
	Knot,
	FeetPerSecondSquared,
	G_Unit,
	MetersPerSecondSquared,
	NewtonsTimeSecond,
	PoundsForceTimesSecond,
	KiloGramsTimesMeterSquared,
	KiloGramsTimesCentiMeterSquared,
	OuncesTimesInchSquared,
	PoundsTimesInchSquared,
	PoundsTimesFeetSquared,
	PoundsForceTimesFeetTimesSecondSquared,
	KiloGramsPerMeterSquared,
	GramsPerCentimeterSquared,
	GramsPerMeterSquared,
	KiloGramsPerMeter,
	OuncesPerInchSquared,
	OuncesPerFeetSquared,
	PoundsPerFeetSquared,
	GramsPerCentimeterCubed,
	KiloGramsPerDeciMeterCubed,
	KiloGramsPerMeterCubed,
	OuncesPerInchCubed,
	OuncesPerFeetCubed,
	PoundsPerFeetCubed,
	MilliMetersSquared,
	CentiMetersSquared,
	InchesSquared,
	FeetSquared,
	MetersSquared,
	Degrees,
	Radians,
	Degrees_Celsius,
	Degrees_Fahrenheit,
	Kelvin,
	Arcminutes,
	DegreesPerSecond,
	RadPerSecond,
	RevolutionsPerSecond,
	RevolutionsPerMinute,
	Bars,
	Pascals,
	Atmospheres,
	MillimeterOfMercury,
	InchesOfMercury,
	PoundsPerSquareInch,  //psi
	Hour,
	Minute,
	Second,
	GramsPerMeter,
	OuncesPerFeet,
	Unknown
};

[[nodiscard]] constexpr f32 to_si(f32 valueue, Type type);
[[nodiscard]] constexpr f32 to(f32 valueue, Type type1, Type type2);

[[nodiscard]] Type type_from_string(char const* string);
[[nodiscard]] char const* type_to_string(Type type);

[[nodiscard]] bool is_unit_of_length(Type unit);
[[nodiscard]] bool is_unit_of_mass(Type unit);
[[nodiscard]] bool is_unit_of_force(Type unit);
[[nodiscard]] bool is_unit_of_velocity(Type unit);
[[nodiscard]] bool is_unit_of_acceleration(Type unit);
[[nodiscard]] bool is_unit_of_impulse(Type unit);
[[nodiscard]] bool is_unit_of_moment_of_inertia(Type unit);
[[nodiscard]] bool is_unit_of_surface_density(Type unit);
[[nodiscard]] bool is_unit_of_bulk_density(Type unit);
[[nodiscard]] bool is_unit_of_area(Type unit);
[[nodiscard]] bool is_unit_of_angle(Type unit);
[[nodiscard]] bool is_unit_of_time(Type unit);
[[nodiscard]] bool is_unit_of_roll_rate(Type unit);
[[nodiscard]] bool is_unit_of_pressure(Type unit);
[[nodiscard]] bool is_unit_of_denisty_line(Type unit);
[[nodiscard]] bool is_unit_of_temperature(Type unit);

#define UNIT_CONVERSION_PROTOTYPE(A, B) \
	[[nodiscard]] constexpr f32 A##_to_##B(f32 value); \
	[[nodiscard]] constexpr f32 B##_to_##A(f32 value);

// SI conversions
UNIT_CONVERSION_PROTOTYPE(Kilo,  Standard)
UNIT_CONVERSION_PROTOTYPE(Hecto, Standard)
UNIT_CONVERSION_PROTOTYPE(Deca,  Standard)
UNIT_CONVERSION_PROTOTYPE(Deci,  Standard)
UNIT_CONVERSION_PROTOTYPE(Centi, Standard)
UNIT_CONVERSION_PROTOTYPE(Milli, Standard)
UNIT_CONVERSION_PROTOTYPE(Mikro, Standard)

// Length
UNIT_CONVERSION_PROTOTYPE(MilliMeters,   Meters)
UNIT_CONVERSION_PROTOTYPE(CentiMeters,   Meters)
UNIT_CONVERSION_PROTOTYPE(KiloMeters,    Meters)
UNIT_CONVERSION_PROTOTYPE(DeciMeters,    Meters)
UNIT_CONVERSION_PROTOTYPE(MikroMeters,   Meters)
UNIT_CONVERSION_PROTOTYPE(InchesOver64,  Meters)
UNIT_CONVERSION_PROTOTYPE(Inches,        Meters)
UNIT_CONVERSION_PROTOTYPE(Feet,          Meters)
UNIT_CONVERSION_PROTOTYPE(Yards,         Meters)
UNIT_CONVERSION_PROTOTYPE(Miles,         Meters)
UNIT_CONVERSION_PROTOTYPE(NauticalMiles, Meters)

// Mass
UNIT_CONVERSION_PROTOTYPE(Grams,  KiloGrams)
UNIT_CONVERSION_PROTOTYPE(Ounces, KiloGrams) //oz
UNIT_CONVERSION_PROTOTYPE(Pounds, KiloGrams) //lb

// Force
UNIT_CONVERSION_PROTOTYPE(PoundsForce,   Newtons)
UNIT_CONVERSION_PROTOTYPE(KilogramForce, Newtons)

// Velocity
UNIT_CONVERSION_PROTOTYPE(KiloMetersPerHour, MetersPerSecond)
UNIT_CONVERSION_PROTOTYPE(FeetPerSecond,     MetersPerSecond)
UNIT_CONVERSION_PROTOTYPE(MilesPerHour,      MetersPerSecond)
UNIT_CONVERSION_PROTOTYPE(Knot,              MetersPerSecond)

// Acceleration
UNIT_CONVERSION_PROTOTYPE(FeetPerSecondSquared, MetersPerSecondSquared)
UNIT_CONVERSION_PROTOTYPE(G_Unit,               MetersPerSecondSquared)

// Impulse
UNIT_CONVERSION_PROTOTYPE(PoundsForceTimesSecond, NewtonsTimeSecond)

// Moment of inertia
UNIT_CONVERSION_PROTOTYPE(KiloGramsTimesCentiMeterSquared,        KiloGramsTimesMeterSquared)
UNIT_CONVERSION_PROTOTYPE(OuncesTimesInchSquared,                 KiloGramsTimesMeterSquared)
UNIT_CONVERSION_PROTOTYPE(PoundsTimesInchSquared,                 KiloGramsTimesMeterSquared)
UNIT_CONVERSION_PROTOTYPE(PoundsTimesFeetSquared,                 KiloGramsTimesMeterSquared)
UNIT_CONVERSION_PROTOTYPE(PoundsForceTimesFeetTimesSecondSquared, KiloGramsTimesMeterSquared)

// Surface density
UNIT_CONVERSION_PROTOTYPE(GramsPerCentimeterSquared, KiloGramsPerMeterSquared)
UNIT_CONVERSION_PROTOTYPE(GramsPerMeterSquared,      KiloGramsPerMeterSquared)
UNIT_CONVERSION_PROTOTYPE(OuncesPerInchSquared,      KiloGramsPerMeterSquared)
UNIT_CONVERSION_PROTOTYPE(OuncesPerFeetSquared,      KiloGramsPerMeterSquared)
UNIT_CONVERSION_PROTOTYPE(PoundsPerFeetSquared,      KiloGramsPerMeterSquared)

// Bulk density
UNIT_CONVERSION_PROTOTYPE(GramsPerCentimeterCubed,    KiloGramsPerMeterCubed)
UNIT_CONVERSION_PROTOTYPE(KiloGramsPerDeciMeterCubed, KiloGramsPerMeterCubed)
UNIT_CONVERSION_PROTOTYPE(OuncesPerInchCubed,         KiloGramsPerMeterCubed)
UNIT_CONVERSION_PROTOTYPE(OuncesPerFeetCubed,         KiloGramsPerMeterCubed)
UNIT_CONVERSION_PROTOTYPE(PoundsPerFeetCubed,         KiloGramsPerMeterCubed)

// Area
UNIT_CONVERSION_PROTOTYPE(MilliMetersSquared, MetersSquared)
UNIT_CONVERSION_PROTOTYPE(CentiMetersSquared, MetersSquared)
UNIT_CONVERSION_PROTOTYPE(InchesSquared,      MetersSquared)
UNIT_CONVERSION_PROTOTYPE(FeetSquared,        MetersSquared)

// Angle
UNIT_CONVERSION_PROTOTYPE(Degrees,    Radians)
UNIT_CONVERSION_PROTOTYPE(Arcminutes, Radians)

//Time 
UNIT_CONVERSION_PROTOTYPE(Hour,   Second)
UNIT_CONVERSION_PROTOTYPE(Minute, Second)

// Roll rate
UNIT_CONVERSION_PROTOTYPE(DegreesPerSecond,     RadPerSecond)
UNIT_CONVERSION_PROTOTYPE(RevolutionsPerSecond, RadPerSecond)
UNIT_CONVERSION_PROTOTYPE(RevolutionsPerMinute, RadPerSecond)

// Pressure
UNIT_CONVERSION_PROTOTYPE(Bars,                Pascals)
UNIT_CONVERSION_PROTOTYPE(Atmospheres,         Pascals)
UNIT_CONVERSION_PROTOTYPE(MillimeterOfMercury, Pascals)
UNIT_CONVERSION_PROTOTYPE(InchesOfMercury,     Pascals)
UNIT_CONVERSION_PROTOTYPE(PoundsPerSquareInch, Pascals) //psi 

//Denisty Line 
UNIT_CONVERSION_PROTOTYPE(GramsPerMeter, KiloGramsPerMeter)
UNIT_CONVERSION_PROTOTYPE(OuncesPerFeet, KiloGramsPerMeter)

#undef UNIT_CONVERSION_PROTOTYPE

//Temperature
[[nodiscard]] constexpr f32 Degrees_Fahrenheit_to_Kelvin(f32 value);
[[nodiscard]] constexpr f32 Degrees_Celsius_to_Kelvin(f32 value);

[[nodiscard]] constexpr f32 Kelvin_to_Degrees_Celsius(f32 value);
[[nodiscard]] constexpr f32 Kelvin_to_Degrees_Fahrenheit(f32 value);

[[nodiscard]] constexpr f32 Degrees_Fahrenheit_to_Degrees_Celsius(f32 value);
[[nodiscard]] constexpr f32 Degrees_Celsius_to_Degrees_Fahrenheit(f32 value);

} // namespace Core::Units
