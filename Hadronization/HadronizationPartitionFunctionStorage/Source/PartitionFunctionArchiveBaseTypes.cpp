// Copyright (c) 2005-2023 Code Synthesis.
//
// This program was generated by CodeSynthesis XSD, an XML Schema to
// C++ data binding compiler.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 2 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
//
// In addition, as a special exception, Code Synthesis gives permission
// to link this program with the Xerces-C++ library (or with modified
// versions of Xerces-C++ that use the same license as Xerces-C++), and
// distribute linked combinations including the two. You must obey the GNU
// General Public License version 2 in all respects for all of the code
// used other than Xerces-C++. If you modify this copy of the program, you
// may extend this exception to your version of the program, but you are
// not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// Furthermore, Code Synthesis makes a special exception for the Free/Libre
// and Open Source Software (FLOSS) which is described in the accompanying
// FLOSSE file.
//

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/pre.hxx>

#include "../Include/PartitionFunctionArchiveBaseTypes.h"

// ChargeConfiguration
//

const ChargeConfiguration::electricCharge_type& ChargeConfiguration::
electricCharge () const
{
  return this->electricCharge_.get ();
}

ChargeConfiguration::electricCharge_type& ChargeConfiguration::
electricCharge ()
{
  return this->electricCharge_.get ();
}

void ChargeConfiguration::
electricCharge (const electricCharge_type& x)
{
  this->electricCharge_.set (x);
}

const ChargeConfiguration::baryonicCharge_type& ChargeConfiguration::
baryonicCharge () const
{
  return this->baryonicCharge_.get ();
}

ChargeConfiguration::baryonicCharge_type& ChargeConfiguration::
baryonicCharge ()
{
  return this->baryonicCharge_.get ();
}

void ChargeConfiguration::
baryonicCharge (const baryonicCharge_type& x)
{
  this->baryonicCharge_.set (x);
}

const ChargeConfiguration::strangeCharge_type& ChargeConfiguration::
strangeCharge () const
{
  return this->strangeCharge_.get ();
}

ChargeConfiguration::strangeCharge_type& ChargeConfiguration::
strangeCharge ()
{
  return this->strangeCharge_.get ();
}

void ChargeConfiguration::
strangeCharge (const strangeCharge_type& x)
{
  this->strangeCharge_.set (x);
}

const ChargeConfiguration::charmCharge_type& ChargeConfiguration::
charmCharge () const
{
  return this->charmCharge_.get ();
}

ChargeConfiguration::charmCharge_type& ChargeConfiguration::
charmCharge ()
{
  return this->charmCharge_.get ();
}

void ChargeConfiguration::
charmCharge (const charmCharge_type& x)
{
  this->charmCharge_.set (x);
}

const ChargeConfiguration::bottomCharge_type& ChargeConfiguration::
bottomCharge () const
{
  return this->bottomCharge_.get ();
}

ChargeConfiguration::bottomCharge_type& ChargeConfiguration::
bottomCharge ()
{
  return this->bottomCharge_.get ();
}

void ChargeConfiguration::
bottomCharge (const bottomCharge_type& x)
{
  this->bottomCharge_.set (x);
}


// MicrocanonicalParameterGridStructure
//

const MicrocanonicalParameterGridStructure::minEnergyDensity_type& MicrocanonicalParameterGridStructure::
minEnergyDensity () const
{
  return this->minEnergyDensity_.get ();
}

MicrocanonicalParameterGridStructure::minEnergyDensity_type& MicrocanonicalParameterGridStructure::
minEnergyDensity ()
{
  return this->minEnergyDensity_.get ();
}

void MicrocanonicalParameterGridStructure::
minEnergyDensity (const minEnergyDensity_type& x)
{
  this->minEnergyDensity_.set (x);
}

const MicrocanonicalParameterGridStructure::maxEnergyDensity_type& MicrocanonicalParameterGridStructure::
maxEnergyDensity () const
{
  return this->maxEnergyDensity_.get ();
}

MicrocanonicalParameterGridStructure::maxEnergyDensity_type& MicrocanonicalParameterGridStructure::
maxEnergyDensity ()
{
  return this->maxEnergyDensity_.get ();
}

void MicrocanonicalParameterGridStructure::
maxEnergyDensity (const maxEnergyDensity_type& x)
{
  this->maxEnergyDensity_.set (x);
}

const MicrocanonicalParameterGridStructure::numberOfEnergyDensityValues_type& MicrocanonicalParameterGridStructure::
numberOfEnergyDensityValues () const
{
  return this->numberOfEnergyDensityValues_.get ();
}

MicrocanonicalParameterGridStructure::numberOfEnergyDensityValues_type& MicrocanonicalParameterGridStructure::
numberOfEnergyDensityValues ()
{
  return this->numberOfEnergyDensityValues_.get ();
}

void MicrocanonicalParameterGridStructure::
numberOfEnergyDensityValues (const numberOfEnergyDensityValues_type& x)
{
  this->numberOfEnergyDensityValues_.set (x);
}

const MicrocanonicalParameterGridStructure::minStrangenessSuppressionParameter_type& MicrocanonicalParameterGridStructure::
minStrangenessSuppressionParameter () const
{
  return this->minStrangenessSuppressionParameter_.get ();
}

MicrocanonicalParameterGridStructure::minStrangenessSuppressionParameter_type& MicrocanonicalParameterGridStructure::
minStrangenessSuppressionParameter ()
{
  return this->minStrangenessSuppressionParameter_.get ();
}

void MicrocanonicalParameterGridStructure::
minStrangenessSuppressionParameter (const minStrangenessSuppressionParameter_type& x)
{
  this->minStrangenessSuppressionParameter_.set (x);
}

const MicrocanonicalParameterGridStructure::maxStrangenessSuppressionParameter_type& MicrocanonicalParameterGridStructure::
maxStrangenessSuppressionParameter () const
{
  return this->maxStrangenessSuppressionParameter_.get ();
}

MicrocanonicalParameterGridStructure::maxStrangenessSuppressionParameter_type& MicrocanonicalParameterGridStructure::
maxStrangenessSuppressionParameter ()
{
  return this->maxStrangenessSuppressionParameter_.get ();
}

void MicrocanonicalParameterGridStructure::
maxStrangenessSuppressionParameter (const maxStrangenessSuppressionParameter_type& x)
{
  this->maxStrangenessSuppressionParameter_.set (x);
}

const MicrocanonicalParameterGridStructure::numberOfStrangenessSuppressionParameterValues_type& MicrocanonicalParameterGridStructure::
numberOfStrangenessSuppressionParameterValues () const
{
  return this->numberOfStrangenessSuppressionParameterValues_.get ();
}

MicrocanonicalParameterGridStructure::numberOfStrangenessSuppressionParameterValues_type& MicrocanonicalParameterGridStructure::
numberOfStrangenessSuppressionParameterValues ()
{
  return this->numberOfStrangenessSuppressionParameterValues_.get ();
}

void MicrocanonicalParameterGridStructure::
numberOfStrangenessSuppressionParameterValues (const numberOfStrangenessSuppressionParameterValues_type& x)
{
  this->numberOfStrangenessSuppressionParameterValues_.set (x);
}


#include <xsd/cxx/xml/dom/parsing-source.hxx>

// ChargeConfiguration
//

ChargeConfiguration::
ChargeConfiguration (const electricCharge_type& electricCharge,
                     const baryonicCharge_type& baryonicCharge,
                     const strangeCharge_type& strangeCharge,
                     const charmCharge_type& charmCharge,
                     const bottomCharge_type& bottomCharge)
: ::xml_schema::type (),
  electricCharge_ (electricCharge, this),
  baryonicCharge_ (baryonicCharge, this),
  strangeCharge_ (strangeCharge, this),
  charmCharge_ (charmCharge, this),
  bottomCharge_ (bottomCharge, this)
{
}

ChargeConfiguration::
ChargeConfiguration (const ChargeConfiguration& x,
                     ::xml_schema::flags f,
                     ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  electricCharge_ (x.electricCharge_, f, this),
  baryonicCharge_ (x.baryonicCharge_, f, this),
  strangeCharge_ (x.strangeCharge_, f, this),
  charmCharge_ (x.charmCharge_, f, this),
  bottomCharge_ (x.bottomCharge_, f, this)
{
}

ChargeConfiguration::
ChargeConfiguration (const ::xercesc::DOMElement& e,
                     ::xml_schema::flags f,
                     ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  electricCharge_ (this),
  baryonicCharge_ (this),
  strangeCharge_ (this),
  charmCharge_ (this),
  bottomCharge_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, false);
    this->parse (p, f);
  }
}

void ChargeConfiguration::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // electricCharge
    //
    if (n.name () == "electricCharge" && n.namespace_ ().empty ())
    {
      if (!electricCharge_.present ())
      {
        this->electricCharge_.set (electricCharge_traits::create (i, f, this));
        continue;
      }
    }

    // baryonicCharge
    //
    if (n.name () == "baryonicCharge" && n.namespace_ ().empty ())
    {
      if (!baryonicCharge_.present ())
      {
        this->baryonicCharge_.set (baryonicCharge_traits::create (i, f, this));
        continue;
      }
    }

    // strangeCharge
    //
    if (n.name () == "strangeCharge" && n.namespace_ ().empty ())
    {
      if (!strangeCharge_.present ())
      {
        this->strangeCharge_.set (strangeCharge_traits::create (i, f, this));
        continue;
      }
    }

    // charmCharge
    //
    if (n.name () == "charmCharge" && n.namespace_ ().empty ())
    {
      if (!charmCharge_.present ())
      {
        this->charmCharge_.set (charmCharge_traits::create (i, f, this));
        continue;
      }
    }

    // bottomCharge
    //
    if (n.name () == "bottomCharge" && n.namespace_ ().empty ())
    {
      if (!bottomCharge_.present ())
      {
        this->bottomCharge_.set (bottomCharge_traits::create (i, f, this));
        continue;
      }
    }

    break;
  }

  if (!electricCharge_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "electricCharge",
      "");
  }

  if (!baryonicCharge_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "baryonicCharge",
      "");
  }

  if (!strangeCharge_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "strangeCharge",
      "");
  }

  if (!charmCharge_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "charmCharge",
      "");
  }

  if (!bottomCharge_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "bottomCharge",
      "");
  }
}

ChargeConfiguration* ChargeConfiguration::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class ChargeConfiguration (*this, f, c);
}

ChargeConfiguration& ChargeConfiguration::
operator= (const ChargeConfiguration& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->electricCharge_ = x.electricCharge_;
    this->baryonicCharge_ = x.baryonicCharge_;
    this->strangeCharge_ = x.strangeCharge_;
    this->charmCharge_ = x.charmCharge_;
    this->bottomCharge_ = x.bottomCharge_;
  }

  return *this;
}

ChargeConfiguration::
~ChargeConfiguration ()
{
}

// MicrocanonicalParameterGridStructure
//

MicrocanonicalParameterGridStructure::
MicrocanonicalParameterGridStructure (const minEnergyDensity_type& minEnergyDensity,
                                      const maxEnergyDensity_type& maxEnergyDensity,
                                      const numberOfEnergyDensityValues_type& numberOfEnergyDensityValues,
                                      const minStrangenessSuppressionParameter_type& minStrangenessSuppressionParameter,
                                      const maxStrangenessSuppressionParameter_type& maxStrangenessSuppressionParameter,
                                      const numberOfStrangenessSuppressionParameterValues_type& numberOfStrangenessSuppressionParameterValues)
: ::xml_schema::type (),
  minEnergyDensity_ (minEnergyDensity, this),
  maxEnergyDensity_ (maxEnergyDensity, this),
  numberOfEnergyDensityValues_ (numberOfEnergyDensityValues, this),
  minStrangenessSuppressionParameter_ (minStrangenessSuppressionParameter, this),
  maxStrangenessSuppressionParameter_ (maxStrangenessSuppressionParameter, this),
  numberOfStrangenessSuppressionParameterValues_ (numberOfStrangenessSuppressionParameterValues, this)
{
}

MicrocanonicalParameterGridStructure::
MicrocanonicalParameterGridStructure (const MicrocanonicalParameterGridStructure& x,
                                      ::xml_schema::flags f,
                                      ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  minEnergyDensity_ (x.minEnergyDensity_, f, this),
  maxEnergyDensity_ (x.maxEnergyDensity_, f, this),
  numberOfEnergyDensityValues_ (x.numberOfEnergyDensityValues_, f, this),
  minStrangenessSuppressionParameter_ (x.minStrangenessSuppressionParameter_, f, this),
  maxStrangenessSuppressionParameter_ (x.maxStrangenessSuppressionParameter_, f, this),
  numberOfStrangenessSuppressionParameterValues_ (x.numberOfStrangenessSuppressionParameterValues_, f, this)
{
}

MicrocanonicalParameterGridStructure::
MicrocanonicalParameterGridStructure (const ::xercesc::DOMElement& e,
                                      ::xml_schema::flags f,
                                      ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  minEnergyDensity_ (this),
  maxEnergyDensity_ (this),
  numberOfEnergyDensityValues_ (this),
  minStrangenessSuppressionParameter_ (this),
  maxStrangenessSuppressionParameter_ (this),
  numberOfStrangenessSuppressionParameterValues_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, false);
    this->parse (p, f);
  }
}

void MicrocanonicalParameterGridStructure::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // minEnergyDensity
    //
    if (n.name () == "minEnergyDensity" && n.namespace_ ().empty ())
    {
      if (!minEnergyDensity_.present ())
      {
        this->minEnergyDensity_.set (minEnergyDensity_traits::create (i, f, this));
        continue;
      }
    }

    // maxEnergyDensity
    //
    if (n.name () == "maxEnergyDensity" && n.namespace_ ().empty ())
    {
      if (!maxEnergyDensity_.present ())
      {
        this->maxEnergyDensity_.set (maxEnergyDensity_traits::create (i, f, this));
        continue;
      }
    }

    // numberOfEnergyDensityValues
    //
    if (n.name () == "numberOfEnergyDensityValues" && n.namespace_ ().empty ())
    {
      if (!numberOfEnergyDensityValues_.present ())
      {
        this->numberOfEnergyDensityValues_.set (numberOfEnergyDensityValues_traits::create (i, f, this));
        continue;
      }
    }

    // minStrangenessSuppressionParameter
    //
    if (n.name () == "minStrangenessSuppressionParameter" && n.namespace_ ().empty ())
    {
      if (!minStrangenessSuppressionParameter_.present ())
      {
        this->minStrangenessSuppressionParameter_.set (minStrangenessSuppressionParameter_traits::create (i, f, this));
        continue;
      }
    }

    // maxStrangenessSuppressionParameter
    //
    if (n.name () == "maxStrangenessSuppressionParameter" && n.namespace_ ().empty ())
    {
      if (!maxStrangenessSuppressionParameter_.present ())
      {
        this->maxStrangenessSuppressionParameter_.set (maxStrangenessSuppressionParameter_traits::create (i, f, this));
        continue;
      }
    }

    // numberOfStrangenessSuppressionParameterValues
    //
    if (n.name () == "numberOfStrangenessSuppressionParameterValues" && n.namespace_ ().empty ())
    {
      if (!numberOfStrangenessSuppressionParameterValues_.present ())
      {
        this->numberOfStrangenessSuppressionParameterValues_.set (numberOfStrangenessSuppressionParameterValues_traits::create (i, f, this));
        continue;
      }
    }

    break;
  }

  if (!minEnergyDensity_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "minEnergyDensity",
      "");
  }

  if (!maxEnergyDensity_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "maxEnergyDensity",
      "");
  }

  if (!numberOfEnergyDensityValues_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "numberOfEnergyDensityValues",
      "");
  }

  if (!minStrangenessSuppressionParameter_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "minStrangenessSuppressionParameter",
      "");
  }

  if (!maxStrangenessSuppressionParameter_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "maxStrangenessSuppressionParameter",
      "");
  }

  if (!numberOfStrangenessSuppressionParameterValues_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "numberOfStrangenessSuppressionParameterValues",
      "");
  }
}

MicrocanonicalParameterGridStructure* MicrocanonicalParameterGridStructure::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class MicrocanonicalParameterGridStructure (*this, f, c);
}

MicrocanonicalParameterGridStructure& MicrocanonicalParameterGridStructure::
operator= (const MicrocanonicalParameterGridStructure& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->minEnergyDensity_ = x.minEnergyDensity_;
    this->maxEnergyDensity_ = x.maxEnergyDensity_;
    this->numberOfEnergyDensityValues_ = x.numberOfEnergyDensityValues_;
    this->minStrangenessSuppressionParameter_ = x.minStrangenessSuppressionParameter_;
    this->maxStrangenessSuppressionParameter_ = x.maxStrangenessSuppressionParameter_;
    this->numberOfStrangenessSuppressionParameterValues_ = x.numberOfStrangenessSuppressionParameterValues_;
  }

  return *this;
}

MicrocanonicalParameterGridStructure::
~MicrocanonicalParameterGridStructure ()
{
}

#include <istream>
#include <xsd/cxx/xml/sax/std-input-source.hxx>
#include <xsd/cxx/tree/error-handler.hxx>

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

