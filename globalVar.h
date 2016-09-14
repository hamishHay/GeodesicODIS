// FILE: globalVar.h
// DESCRIPTION: header file for the GlobalVar and base IGlobalVar class. GlobalVar
//							is a template class, allowing for its main variable to be
//							initialised of any type. Here we use int, double, string, and
// 							bool. It inherits attributes from the base class IGlobalVar, mainly
//							its string ID (used to identify the variable) and type (as a string
//							in order to check the incoming input type).
//
// 	Version 			 			Date 									Programmer
// 		0.1			|				13/09/2016				|				H. Hay				|
// ---- Initial version of ODIS, described in Hay and Matsuyama (2016)

#ifndef GLOBALVAR_H
#define GLOBALVAR_H

#include <string>
#include <iostream>

//Base class
class IGlobalVar
{
public:

	// virtual objects are defined here to be able to identify attributes about each
	// global variable through its base class.
	virtual std::string StringID() = 0;
	virtual bool IsType(std::string) = 0;

	void Added(bool boolValue)
	{
		added = boolValue;
	}

	bool IsAdded()
	{
		return added;
	}

	void SetStringID(std::string n)
	{
		stringID = n;
	};

protected:
	// define variables common to all global variables
	std::string stringID; // variable ID tag
	std::string type; // variable type as a string
	bool added = false;
};

// template class, based on IGlobalVar. The class T entered into the template
// allows us to define the actual type of the global variable stored by a GlobalVar
// instance.
template <class T>
class GlobalVar : public IGlobalVar
{
private:
	void SetType()
	{
		if (std::is_same<T, int>::value) type = "int";
		else if (std::is_same<T, double>::value) type = "double";
		else if (std::is_same<T, bool>::value) type = "bool";
		else if (std::is_same<T, std::string>::value) type = "string";
	}

public:
	GlobalVar()
	{
		SetType(); //Constructor calls function to define class type;
	}

	// virtual attributes inherited from the base class.
	virtual std::string StringID()
	{
		return stringID;
	}

	virtual bool IsType(std::string typeCheck)
	{
		if (type == typeCheck) return true;
		else return false;
	};

	// set value of variable of unknown type T
	void SetValue(T n)
	{
		val = n;
	};

	// return value of variable of unknown type T
	inline T Value(void)
	{
		return val;
	};

protected:
	T val; // define variable of uknown type T
};

#endif
