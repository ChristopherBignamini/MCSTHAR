// Copyright (C) 2005-2010 Code Synthesis Tools CC
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
// In addition, as a special exception, Code Synthesis Tools CC gives
// permission to link this program with the Xerces-C++ library (or with
// modified versions of Xerces-C++ that use the same license as Xerces-C++),
// and distribute linked combinations including the two. You must obey
// the GNU General Public License version 2 in all respects for all of
// the code used other than Xerces-C++. If you modify this copy of the
// program, you may extend this exception to your version of the program,
// but you are not obligated to do so. If you do not wish to do so, delete
// this exception statement from your version.
//
// Furthermore, Code Synthesis Tools CC makes a special exception for
// the Free/Libre and Open Source Software (FLOSS) which is described
// in the accompanying FLOSSE file.
//

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/pre.hxx>

#include "../Include/MCSTHARSetup.h"

namespace MCSTHARSetup
{
  // MCSTHARSetup
  // 

  const MCSTHARSetup::partitionFunctionDataSetPath_type& MCSTHARSetup::
  partitionFunctionDataSetPath () const
  {
    return this->partitionFunctionDataSetPath_.get ();
  }

  MCSTHARSetup::partitionFunctionDataSetPath_type& MCSTHARSetup::
  partitionFunctionDataSetPath ()
  {
    return this->partitionFunctionDataSetPath_.get ();
  }

  void MCSTHARSetup::
  partitionFunctionDataSetPath (const partitionFunctionDataSetPath_type& x)
  {
    this->partitionFunctionDataSetPath_.set (x);
  }

  void MCSTHARSetup::
  partitionFunctionDataSetPath (::std::auto_ptr< partitionFunctionDataSetPath_type > x)
  {
    this->partitionFunctionDataSetPath_.set (x);
  }

  const MCSTHARSetup::hadronDataFileName_type& MCSTHARSetup::
  hadronDataFileName () const
  {
    return this->hadronDataFileName_.get ();
  }

  MCSTHARSetup::hadronDataFileName_type& MCSTHARSetup::
  hadronDataFileName ()
  {
    return this->hadronDataFileName_.get ();
  }

  void MCSTHARSetup::
  hadronDataFileName (const hadronDataFileName_type& x)
  {
    this->hadronDataFileName_.set (x);
  }

  void MCSTHARSetup::
  hadronDataFileName (::std::auto_ptr< hadronDataFileName_type > x)
  {
    this->hadronDataFileName_.set (x);
  }

  const MCSTHARSetup::clusterEnergyDensity_type& MCSTHARSetup::
  clusterEnergyDensity () const
  {
    return this->clusterEnergyDensity_.get ();
  }

  MCSTHARSetup::clusterEnergyDensity_type& MCSTHARSetup::
  clusterEnergyDensity ()
  {
    return this->clusterEnergyDensity_.get ();
  }

  void MCSTHARSetup::
  clusterEnergyDensity (const clusterEnergyDensity_type& x)
  {
    this->clusterEnergyDensity_.set (x);
  }

  const MCSTHARSetup::strangenessSuppressionParameter_type& MCSTHARSetup::
  strangenessSuppressionParameter () const
  {
    return this->strangenessSuppressionParameter_.get ();
  }

  MCSTHARSetup::strangenessSuppressionParameter_type& MCSTHARSetup::
  strangenessSuppressionParameter ()
  {
    return this->strangenessSuppressionParameter_.get ();
  }

  void MCSTHARSetup::
  strangenessSuppressionParameter (const strangenessSuppressionParameter_type& x)
  {
    this->strangenessSuppressionParameter_.set (x);
  }

  const MCSTHARSetup::clusterMergingMinimumMass_optional& MCSTHARSetup::
  clusterMergingMinimumMass () const
  {
    return this->clusterMergingMinimumMass_;
  }

  MCSTHARSetup::clusterMergingMinimumMass_optional& MCSTHARSetup::
  clusterMergingMinimumMass ()
  {
    return this->clusterMergingMinimumMass_;
  }

  void MCSTHARSetup::
  clusterMergingMinimumMass (const clusterMergingMinimumMass_type& x)
  {
    this->clusterMergingMinimumMass_.set (x);
  }

  void MCSTHARSetup::
  clusterMergingMinimumMass (const clusterMergingMinimumMass_optional& x)
  {
    this->clusterMergingMinimumMass_ = x;
  }

  const MCSTHARSetup::charmClusterMergingMinimumMass_optional& MCSTHARSetup::
  charmClusterMergingMinimumMass () const
  {
    return this->charmClusterMergingMinimumMass_;
  }

  MCSTHARSetup::charmClusterMergingMinimumMass_optional& MCSTHARSetup::
  charmClusterMergingMinimumMass ()
  {
    return this->charmClusterMergingMinimumMass_;
  }

  void MCSTHARSetup::
  charmClusterMergingMinimumMass (const charmClusterMergingMinimumMass_type& x)
  {
    this->charmClusterMergingMinimumMass_.set (x);
  }

  void MCSTHARSetup::
  charmClusterMergingMinimumMass (const charmClusterMergingMinimumMass_optional& x)
  {
    this->charmClusterMergingMinimumMass_ = x;
  }

  const MCSTHARSetup::bottomClusterMergingMinimumMass_optional& MCSTHARSetup::
  bottomClusterMergingMinimumMass () const
  {
    return this->bottomClusterMergingMinimumMass_;
  }

  MCSTHARSetup::bottomClusterMergingMinimumMass_optional& MCSTHARSetup::
  bottomClusterMergingMinimumMass ()
  {
    return this->bottomClusterMergingMinimumMass_;
  }

  void MCSTHARSetup::
  bottomClusterMergingMinimumMass (const bottomClusterMergingMinimumMass_type& x)
  {
    this->bottomClusterMergingMinimumMass_.set (x);
  }

  void MCSTHARSetup::
  bottomClusterMergingMinimumMass (const bottomClusterMergingMinimumMass_optional& x)
  {
    this->bottomClusterMergingMinimumMass_ = x;
  }

  const MCSTHARSetup::randomNumberGeneratorSeed_optional& MCSTHARSetup::
  randomNumberGeneratorSeed () const
  {
    return this->randomNumberGeneratorSeed_;
  }

  MCSTHARSetup::randomNumberGeneratorSeed_optional& MCSTHARSetup::
  randomNumberGeneratorSeed ()
  {
    return this->randomNumberGeneratorSeed_;
  }

  void MCSTHARSetup::
  randomNumberGeneratorSeed (const randomNumberGeneratorSeed_type& x)
  {
    this->randomNumberGeneratorSeed_.set (x);
  }

  void MCSTHARSetup::
  randomNumberGeneratorSeed (const randomNumberGeneratorSeed_optional& x)
  {
    this->randomNumberGeneratorSeed_ = x;
  }
}

#include <xsd/cxx/xml/dom/parsing-source.hxx>

namespace MCSTHARSetup
{
  // MCSTHARSetup
  //

  MCSTHARSetup::
  MCSTHARSetup (const partitionFunctionDataSetPath_type& partitionFunctionDataSetPath,
                const hadronDataFileName_type& hadronDataFileName,
                const clusterEnergyDensity_type& clusterEnergyDensity,
                const strangenessSuppressionParameter_type& strangenessSuppressionParameter)
  : ::xml_schema::type (),
    partitionFunctionDataSetPath_ (partitionFunctionDataSetPath, ::xml_schema::flags (), this),
    hadronDataFileName_ (hadronDataFileName, ::xml_schema::flags (), this),
    clusterEnergyDensity_ (clusterEnergyDensity, ::xml_schema::flags (), this),
    strangenessSuppressionParameter_ (strangenessSuppressionParameter, ::xml_schema::flags (), this),
    clusterMergingMinimumMass_ (::xml_schema::flags (), this),
    charmClusterMergingMinimumMass_ (::xml_schema::flags (), this),
    bottomClusterMergingMinimumMass_ (::xml_schema::flags (), this),
    randomNumberGeneratorSeed_ (::xml_schema::flags (), this)
  {
  }

  MCSTHARSetup::
  MCSTHARSetup (const MCSTHARSetup& x,
                ::xml_schema::flags f,
                ::xml_schema::container* c)
  : ::xml_schema::type (x, f, c),
    partitionFunctionDataSetPath_ (x.partitionFunctionDataSetPath_, f, this),
    hadronDataFileName_ (x.hadronDataFileName_, f, this),
    clusterEnergyDensity_ (x.clusterEnergyDensity_, f, this),
    strangenessSuppressionParameter_ (x.strangenessSuppressionParameter_, f, this),
    clusterMergingMinimumMass_ (x.clusterMergingMinimumMass_, f, this),
    charmClusterMergingMinimumMass_ (x.charmClusterMergingMinimumMass_, f, this),
    bottomClusterMergingMinimumMass_ (x.bottomClusterMergingMinimumMass_, f, this),
    randomNumberGeneratorSeed_ (x.randomNumberGeneratorSeed_, f, this)
  {
  }

  MCSTHARSetup::
  MCSTHARSetup (const ::xercesc::DOMElement& e,
                ::xml_schema::flags f,
                ::xml_schema::container* c)
  : ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
    partitionFunctionDataSetPath_ (f, this),
    hadronDataFileName_ (f, this),
    clusterEnergyDensity_ (f, this),
    strangenessSuppressionParameter_ (f, this),
    clusterMergingMinimumMass_ (f, this),
    charmClusterMergingMinimumMass_ (f, this),
    bottomClusterMergingMinimumMass_ (f, this),
    randomNumberGeneratorSeed_ (f, this)
  {
    if ((f & ::xml_schema::flags::base) == 0)
    {
      ::xsd::cxx::xml::dom::parser< char > p (e, true, false);
      this->parse (p, f);
    }
  }

  void MCSTHARSetup::
  parse (::xsd::cxx::xml::dom::parser< char >& p,
         ::xml_schema::flags f)
  {
    for (; p.more_elements (); p.next_element ())
    {
      const ::xercesc::DOMElement& i (p.cur_element ());
      const ::xsd::cxx::xml::qualified_name< char > n (
        ::xsd::cxx::xml::dom::name< char > (i));

      // partitionFunctionDataSetPath
      //
      if (n.name () == "partitionFunctionDataSetPath" && n.namespace_ ().empty ())
      {
        ::std::auto_ptr< partitionFunctionDataSetPath_type > r (
          partitionFunctionDataSetPath_traits::create (i, f, this));

        if (!partitionFunctionDataSetPath_.present ())
        {
          this->partitionFunctionDataSetPath_.set (r);
          continue;
        }
      }

      // hadronDataFileName
      //
      if (n.name () == "hadronDataFileName" && n.namespace_ ().empty ())
      {
        ::std::auto_ptr< hadronDataFileName_type > r (
          hadronDataFileName_traits::create (i, f, this));

        if (!hadronDataFileName_.present ())
        {
          this->hadronDataFileName_.set (r);
          continue;
        }
      }

      // clusterEnergyDensity
      //
      if (n.name () == "clusterEnergyDensity" && n.namespace_ ().empty ())
      {
        if (!clusterEnergyDensity_.present ())
        {
          this->clusterEnergyDensity_.set (clusterEnergyDensity_traits::create (i, f, this));
          continue;
        }
      }

      // strangenessSuppressionParameter
      //
      if (n.name () == "strangenessSuppressionParameter" && n.namespace_ ().empty ())
      {
        if (!strangenessSuppressionParameter_.present ())
        {
          this->strangenessSuppressionParameter_.set (strangenessSuppressionParameter_traits::create (i, f, this));
          continue;
        }
      }

      // clusterMergingMinimumMass
      //
      if (n.name () == "clusterMergingMinimumMass" && n.namespace_ ().empty ())
      {
        if (!this->clusterMergingMinimumMass_)
        {
          this->clusterMergingMinimumMass_.set (clusterMergingMinimumMass_traits::create (i, f, this));
          continue;
        }
      }

      // charmClusterMergingMinimumMass
      //
      if (n.name () == "charmClusterMergingMinimumMass" && n.namespace_ ().empty ())
      {
        if (!this->charmClusterMergingMinimumMass_)
        {
          this->charmClusterMergingMinimumMass_.set (charmClusterMergingMinimumMass_traits::create (i, f, this));
          continue;
        }
      }

      // bottomClusterMergingMinimumMass
      //
      if (n.name () == "bottomClusterMergingMinimumMass" && n.namespace_ ().empty ())
      {
        if (!this->bottomClusterMergingMinimumMass_)
        {
          this->bottomClusterMergingMinimumMass_.set (bottomClusterMergingMinimumMass_traits::create (i, f, this));
          continue;
        }
      }

      // randomNumberGeneratorSeed
      //
      if (n.name () == "randomNumberGeneratorSeed" && n.namespace_ ().empty ())
      {
        if (!this->randomNumberGeneratorSeed_)
        {
          this->randomNumberGeneratorSeed_.set (randomNumberGeneratorSeed_traits::create (i, f, this));
          continue;
        }
      }

      break;
    }

    if (!partitionFunctionDataSetPath_.present ())
    {
      throw ::xsd::cxx::tree::expected_element< char > (
        "partitionFunctionDataSetPath",
        "");
    }

    if (!hadronDataFileName_.present ())
    {
      throw ::xsd::cxx::tree::expected_element< char > (
        "hadronDataFileName",
        "");
    }

    if (!clusterEnergyDensity_.present ())
    {
      throw ::xsd::cxx::tree::expected_element< char > (
        "clusterEnergyDensity",
        "");
    }

    if (!strangenessSuppressionParameter_.present ())
    {
      throw ::xsd::cxx::tree::expected_element< char > (
        "strangenessSuppressionParameter",
        "");
    }
  }

  MCSTHARSetup* MCSTHARSetup::
  _clone (::xml_schema::flags f,
          ::xml_schema::container* c) const
  {
    return new class MCSTHARSetup (*this, f, c);
  }

  MCSTHARSetup::
  ~MCSTHARSetup ()
  {
  }
}

#include <istream>
#include <xsd/cxx/xml/sax/std-input-source.hxx>
#include <xsd/cxx/tree/error-handler.hxx>

namespace MCSTHARSetup
{
  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (const ::std::string& u,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::auto_initializer i (
      (f & ::xml_schema::flags::dont_initialize) == 0,
      (f & ::xml_schema::flags::keep_dom) == 0);

    ::xsd::cxx::tree::error_handler< char > h;

    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      ::xsd::cxx::xml::dom::parse< char > (
        u, h, p, f));

    h.throw_if_failed< ::xsd::cxx::tree::parsing< char > > ();

    ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup > r (
      ::MCSTHARSetup::MCSTHARSetup_ (
        d, f | ::xml_schema::flags::own_dom, p));

    return r;
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (const ::std::string& u,
                 ::xml_schema::error_handler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::auto_initializer i (
      (f & ::xml_schema::flags::dont_initialize) == 0,
      (f & ::xml_schema::flags::keep_dom) == 0);

    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      ::xsd::cxx::xml::dom::parse< char > (
        u, h, p, f));

    if (!d.get ())
      throw ::xsd::cxx::tree::parsing< char > ();

    ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup > r (
      ::MCSTHARSetup::MCSTHARSetup_ (
        d, f | ::xml_schema::flags::own_dom, p));

    return r;
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (const ::std::string& u,
                 ::xercesc::DOMErrorHandler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      ::xsd::cxx::xml::dom::parse< char > (
        u, h, p, f));

    if (!d.get ())
      throw ::xsd::cxx::tree::parsing< char > ();

    ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup > r (
      ::MCSTHARSetup::MCSTHARSetup_ (
        d, f | ::xml_schema::flags::own_dom, p));

    return r;
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (::std::istream& is,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::auto_initializer i (
      (f & ::xml_schema::flags::dont_initialize) == 0,
      (f & ::xml_schema::flags::keep_dom) == 0);

    ::xsd::cxx::xml::sax::std_input_source isrc (is);
    return ::MCSTHARSetup::MCSTHARSetup_ (isrc, f, p);
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (::std::istream& is,
                 ::xml_schema::error_handler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::auto_initializer i (
      (f & ::xml_schema::flags::dont_initialize) == 0,
      (f & ::xml_schema::flags::keep_dom) == 0);

    ::xsd::cxx::xml::sax::std_input_source isrc (is);
    return ::MCSTHARSetup::MCSTHARSetup_ (isrc, h, f, p);
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (::std::istream& is,
                 ::xercesc::DOMErrorHandler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::sax::std_input_source isrc (is);
    return ::MCSTHARSetup::MCSTHARSetup_ (isrc, h, f, p);
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (::std::istream& is,
                 const ::std::string& sid,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::auto_initializer i (
      (f & ::xml_schema::flags::dont_initialize) == 0,
      (f & ::xml_schema::flags::keep_dom) == 0);

    ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
    return ::MCSTHARSetup::MCSTHARSetup_ (isrc, f, p);
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (::std::istream& is,
                 const ::std::string& sid,
                 ::xml_schema::error_handler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::auto_initializer i (
      (f & ::xml_schema::flags::dont_initialize) == 0,
      (f & ::xml_schema::flags::keep_dom) == 0);

    ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
    return ::MCSTHARSetup::MCSTHARSetup_ (isrc, h, f, p);
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (::std::istream& is,
                 const ::std::string& sid,
                 ::xercesc::DOMErrorHandler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
    return ::MCSTHARSetup::MCSTHARSetup_ (isrc, h, f, p);
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (::xercesc::InputSource& i,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    ::xsd::cxx::tree::error_handler< char > h;

    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      ::xsd::cxx::xml::dom::parse< char > (
        i, h, p, f));

    h.throw_if_failed< ::xsd::cxx::tree::parsing< char > > ();

    ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup > r (
      ::MCSTHARSetup::MCSTHARSetup_ (
        d, f | ::xml_schema::flags::own_dom, p));

    return r;
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (::xercesc::InputSource& i,
                 ::xml_schema::error_handler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      ::xsd::cxx::xml::dom::parse< char > (
        i, h, p, f));

    if (!d.get ())
      throw ::xsd::cxx::tree::parsing< char > ();

    ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup > r (
      ::MCSTHARSetup::MCSTHARSetup_ (
        d, f | ::xml_schema::flags::own_dom, p));

    return r;
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (::xercesc::InputSource& i,
                 ::xercesc::DOMErrorHandler& h,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > d (
      ::xsd::cxx::xml::dom::parse< char > (
        i, h, p, f));

    if (!d.get ())
      throw ::xsd::cxx::tree::parsing< char > ();

    ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup > r (
      ::MCSTHARSetup::MCSTHARSetup_ (
        d, f | ::xml_schema::flags::own_dom, p));

    return r;
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (const ::xercesc::DOMDocument& d,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties& p)
  {
    if (f & ::xml_schema::flags::keep_dom)
    {
      ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > c (
        static_cast< ::xercesc::DOMDocument* > (d.cloneNode (true)));

      ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup > r (
        ::MCSTHARSetup::MCSTHARSetup_ (
          c, f | ::xml_schema::flags::own_dom, p));

      return r;
    }

    const ::xercesc::DOMElement& e (*d.getDocumentElement ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (e));

    if (n.name () == "MCSTHARSetup" &&
        n.namespace_ () == "")
    {
      ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup > r (
        ::xsd::cxx::tree::traits< ::MCSTHARSetup::MCSTHARSetup, char >::create (
          e, f, 0));
      return r;
    }

    throw ::xsd::cxx::tree::unexpected_element < char > (
      n.name (),
      n.namespace_ (),
      "MCSTHARSetup",
      "");
  }

  ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup >
  MCSTHARSetup_ (::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument >& d,
                 ::xml_schema::flags f,
                 const ::xml_schema::properties&)
  {
    ::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument > c (
      ((f & ::xml_schema::flags::keep_dom) &&
       !(f & ::xml_schema::flags::own_dom))
      ? static_cast< ::xercesc::DOMDocument* > (d->cloneNode (true))
      : 0);

    ::xercesc::DOMDocument& doc (c.get () ? *c : *d);
    const ::xercesc::DOMElement& e (*doc.getDocumentElement ());

    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (e));

    if (f & ::xml_schema::flags::keep_dom)
      doc.setUserData (::xml_schema::dom::tree_node_key,
                       (c.get () ? &c : &d),
                       0);

    if (n.name () == "MCSTHARSetup" &&
        n.namespace_ () == "")
    {
      ::std::auto_ptr< ::MCSTHARSetup::MCSTHARSetup > r (
        ::xsd::cxx::tree::traits< ::MCSTHARSetup::MCSTHARSetup, char >::create (
          e, f, 0));
      return r;
    }

    throw ::xsd::cxx::tree::unexpected_element < char > (
      n.name (),
      n.namespace_ (),
      "MCSTHARSetup",
      "");
  }
}

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

