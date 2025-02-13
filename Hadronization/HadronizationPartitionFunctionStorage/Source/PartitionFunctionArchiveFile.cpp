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

#include "../Include/PartitionFunctionArchiveFile.h"

// SingleChargeConfigurationArchive
//

const SingleChargeConfigurationArchive::chargeConfiguration_type& SingleChargeConfigurationArchive::
chargeConfiguration () const
{
  return this->chargeConfiguration_.get ();
}

SingleChargeConfigurationArchive::chargeConfiguration_type& SingleChargeConfigurationArchive::
chargeConfiguration ()
{
  return this->chargeConfiguration_.get ();
}

void SingleChargeConfigurationArchive::
chargeConfiguration (const chargeConfiguration_type& x)
{
  this->chargeConfiguration_.set (x);
}

void SingleChargeConfigurationArchive::
chargeConfiguration (::std::unique_ptr< chargeConfiguration_type > x)
{
  this->chargeConfiguration_.set (std::move (x));
}

const SingleChargeConfigurationArchive::folderName_type& SingleChargeConfigurationArchive::
folderName () const
{
  return this->folderName_.get ();
}

SingleChargeConfigurationArchive::folderName_type& SingleChargeConfigurationArchive::
folderName ()
{
  return this->folderName_.get ();
}

void SingleChargeConfigurationArchive::
folderName (const folderName_type& x)
{
  this->folderName_.set (x);
}

void SingleChargeConfigurationArchive::
folderName (::std::unique_ptr< folderName_type > x)
{
  this->folderName_.set (std::move (x));
}


// SingleChargeConfigurationArchiveList
//

const SingleChargeConfigurationArchiveList::singleChargeConfigurationArchive_sequence& SingleChargeConfigurationArchiveList::
singleChargeConfigurationArchive () const
{
  return this->singleChargeConfigurationArchive_;
}

SingleChargeConfigurationArchiveList::singleChargeConfigurationArchive_sequence& SingleChargeConfigurationArchiveList::
singleChargeConfigurationArchive ()
{
  return this->singleChargeConfigurationArchive_;
}

void SingleChargeConfigurationArchiveList::
singleChargeConfigurationArchive (const singleChargeConfigurationArchive_sequence& s)
{
  this->singleChargeConfigurationArchive_ = s;
}

const SingleChargeConfigurationArchiveList::length_type& SingleChargeConfigurationArchiveList::
length () const
{
  return this->length_.get ();
}

SingleChargeConfigurationArchiveList::length_type& SingleChargeConfigurationArchiveList::
length ()
{
  return this->length_.get ();
}

void SingleChargeConfigurationArchiveList::
length (const length_type& x)
{
  this->length_.set (x);
}


// PartitionFunctionArchiveFile
//

const PartitionFunctionArchiveFile::microcanonicalParameterGridStructure_type& PartitionFunctionArchiveFile::
microcanonicalParameterGridStructure () const
{
  return this->microcanonicalParameterGridStructure_.get ();
}

PartitionFunctionArchiveFile::microcanonicalParameterGridStructure_type& PartitionFunctionArchiveFile::
microcanonicalParameterGridStructure ()
{
  return this->microcanonicalParameterGridStructure_.get ();
}

void PartitionFunctionArchiveFile::
microcanonicalParameterGridStructure (const microcanonicalParameterGridStructure_type& x)
{
  this->microcanonicalParameterGridStructure_.set (x);
}

void PartitionFunctionArchiveFile::
microcanonicalParameterGridStructure (::std::unique_ptr< microcanonicalParameterGridStructure_type > x)
{
  this->microcanonicalParameterGridStructure_.set (std::move (x));
}

const PartitionFunctionArchiveFile::singleChargeConfigurationArchiveList_type& PartitionFunctionArchiveFile::
singleChargeConfigurationArchiveList () const
{
  return this->singleChargeConfigurationArchiveList_.get ();
}

PartitionFunctionArchiveFile::singleChargeConfigurationArchiveList_type& PartitionFunctionArchiveFile::
singleChargeConfigurationArchiveList ()
{
  return this->singleChargeConfigurationArchiveList_.get ();
}

void PartitionFunctionArchiveFile::
singleChargeConfigurationArchiveList (const singleChargeConfigurationArchiveList_type& x)
{
  this->singleChargeConfigurationArchiveList_.set (x);
}

void PartitionFunctionArchiveFile::
singleChargeConfigurationArchiveList (::std::unique_ptr< singleChargeConfigurationArchiveList_type > x)
{
  this->singleChargeConfigurationArchiveList_.set (std::move (x));
}


#include <xsd/cxx/xml/dom/parsing-source.hxx>

// SingleChargeConfigurationArchive
//

SingleChargeConfigurationArchive::
SingleChargeConfigurationArchive (const chargeConfiguration_type& chargeConfiguration,
                                  const folderName_type& folderName)
: ::xml_schema::type (),
  chargeConfiguration_ (chargeConfiguration, this),
  folderName_ (folderName, this)
{
}

SingleChargeConfigurationArchive::
SingleChargeConfigurationArchive (::std::unique_ptr< chargeConfiguration_type > chargeConfiguration,
                                  const folderName_type& folderName)
: ::xml_schema::type (),
  chargeConfiguration_ (std::move (chargeConfiguration), this),
  folderName_ (folderName, this)
{
}

SingleChargeConfigurationArchive::
SingleChargeConfigurationArchive (const SingleChargeConfigurationArchive& x,
                                  ::xml_schema::flags f,
                                  ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  chargeConfiguration_ (x.chargeConfiguration_, f, this),
  folderName_ (x.folderName_, f, this)
{
}

SingleChargeConfigurationArchive::
SingleChargeConfigurationArchive (const ::xercesc::DOMElement& e,
                                  ::xml_schema::flags f,
                                  ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  chargeConfiguration_ (this),
  folderName_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, false);
    this->parse (p, f);
  }
}

void SingleChargeConfigurationArchive::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // chargeConfiguration
    //
    if (n.name () == "chargeConfiguration" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< chargeConfiguration_type > r (
        chargeConfiguration_traits::create (i, f, this));

      if (!chargeConfiguration_.present ())
      {
        this->chargeConfiguration_.set (::std::move (r));
        continue;
      }
    }

    // folderName
    //
    if (n.name () == "folderName" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< folderName_type > r (
        folderName_traits::create (i, f, this));

      if (!folderName_.present ())
      {
        this->folderName_.set (::std::move (r));
        continue;
      }
    }

    break;
  }

  if (!chargeConfiguration_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "chargeConfiguration",
      "");
  }

  if (!folderName_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "folderName",
      "");
  }
}

SingleChargeConfigurationArchive* SingleChargeConfigurationArchive::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class SingleChargeConfigurationArchive (*this, f, c);
}

SingleChargeConfigurationArchive& SingleChargeConfigurationArchive::
operator= (const SingleChargeConfigurationArchive& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->chargeConfiguration_ = x.chargeConfiguration_;
    this->folderName_ = x.folderName_;
  }

  return *this;
}

SingleChargeConfigurationArchive::
~SingleChargeConfigurationArchive ()
{
}

// SingleChargeConfigurationArchiveList
//

SingleChargeConfigurationArchiveList::
SingleChargeConfigurationArchiveList (const length_type& length)
: ::xml_schema::type (),
  singleChargeConfigurationArchive_ (this),
  length_ (length, this)
{
}

SingleChargeConfigurationArchiveList::
SingleChargeConfigurationArchiveList (const SingleChargeConfigurationArchiveList& x,
                                      ::xml_schema::flags f,
                                      ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  singleChargeConfigurationArchive_ (x.singleChargeConfigurationArchive_, f, this),
  length_ (x.length_, f, this)
{
}

SingleChargeConfigurationArchiveList::
SingleChargeConfigurationArchiveList (const ::xercesc::DOMElement& e,
                                      ::xml_schema::flags f,
                                      ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  singleChargeConfigurationArchive_ (this),
  length_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, true);
    this->parse (p, f);
  }
}

void SingleChargeConfigurationArchiveList::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // singleChargeConfigurationArchive
    //
    if (n.name () == "singleChargeConfigurationArchive" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< singleChargeConfigurationArchive_type > r (
        singleChargeConfigurationArchive_traits::create (i, f, this));

      this->singleChargeConfigurationArchive_.push_back (::std::move (r));
      continue;
    }

    break;
  }

  while (p.more_attributes ())
  {
    const ::xercesc::DOMAttr& i (p.next_attribute ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    if (n.name () == "length" && n.namespace_ ().empty ())
    {
      this->length_.set (length_traits::create (i, f, this));
      continue;
    }
  }

  if (!length_.present ())
  {
    throw ::xsd::cxx::tree::expected_attribute< char > (
      "length",
      "");
  }
}

SingleChargeConfigurationArchiveList* SingleChargeConfigurationArchiveList::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class SingleChargeConfigurationArchiveList (*this, f, c);
}

SingleChargeConfigurationArchiveList& SingleChargeConfigurationArchiveList::
operator= (const SingleChargeConfigurationArchiveList& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->singleChargeConfigurationArchive_ = x.singleChargeConfigurationArchive_;
    this->length_ = x.length_;
  }

  return *this;
}

SingleChargeConfigurationArchiveList::
~SingleChargeConfigurationArchiveList ()
{
}

// PartitionFunctionArchiveFile
//

PartitionFunctionArchiveFile::
PartitionFunctionArchiveFile (const microcanonicalParameterGridStructure_type& microcanonicalParameterGridStructure,
                              const singleChargeConfigurationArchiveList_type& singleChargeConfigurationArchiveList)
: ::xml_schema::type (),
  microcanonicalParameterGridStructure_ (microcanonicalParameterGridStructure, this),
  singleChargeConfigurationArchiveList_ (singleChargeConfigurationArchiveList, this)
{
}

PartitionFunctionArchiveFile::
PartitionFunctionArchiveFile (::std::unique_ptr< microcanonicalParameterGridStructure_type > microcanonicalParameterGridStructure,
                              ::std::unique_ptr< singleChargeConfigurationArchiveList_type > singleChargeConfigurationArchiveList)
: ::xml_schema::type (),
  microcanonicalParameterGridStructure_ (std::move (microcanonicalParameterGridStructure), this),
  singleChargeConfigurationArchiveList_ (std::move (singleChargeConfigurationArchiveList), this)
{
}

PartitionFunctionArchiveFile::
PartitionFunctionArchiveFile (const PartitionFunctionArchiveFile& x,
                              ::xml_schema::flags f,
                              ::xml_schema::container* c)
: ::xml_schema::type (x, f, c),
  microcanonicalParameterGridStructure_ (x.microcanonicalParameterGridStructure_, f, this),
  singleChargeConfigurationArchiveList_ (x.singleChargeConfigurationArchiveList_, f, this)
{
}

PartitionFunctionArchiveFile::
PartitionFunctionArchiveFile (const ::xercesc::DOMElement& e,
                              ::xml_schema::flags f,
                              ::xml_schema::container* c)
: ::xml_schema::type (e, f | ::xml_schema::flags::base, c),
  microcanonicalParameterGridStructure_ (this),
  singleChargeConfigurationArchiveList_ (this)
{
  if ((f & ::xml_schema::flags::base) == 0)
  {
    ::xsd::cxx::xml::dom::parser< char > p (e, true, false, false);
    this->parse (p, f);
  }
}

void PartitionFunctionArchiveFile::
parse (::xsd::cxx::xml::dom::parser< char >& p,
       ::xml_schema::flags f)
{
  for (; p.more_content (); p.next_content (false))
  {
    const ::xercesc::DOMElement& i (p.cur_element ());
    const ::xsd::cxx::xml::qualified_name< char > n (
      ::xsd::cxx::xml::dom::name< char > (i));

    // microcanonicalParameterGridStructure
    //
    if (n.name () == "microcanonicalParameterGridStructure" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< microcanonicalParameterGridStructure_type > r (
        microcanonicalParameterGridStructure_traits::create (i, f, this));

      if (!microcanonicalParameterGridStructure_.present ())
      {
        this->microcanonicalParameterGridStructure_.set (::std::move (r));
        continue;
      }
    }

    // singleChargeConfigurationArchiveList
    //
    if (n.name () == "singleChargeConfigurationArchiveList" && n.namespace_ ().empty ())
    {
      ::std::unique_ptr< singleChargeConfigurationArchiveList_type > r (
        singleChargeConfigurationArchiveList_traits::create (i, f, this));

      if (!singleChargeConfigurationArchiveList_.present ())
      {
        this->singleChargeConfigurationArchiveList_.set (::std::move (r));
        continue;
      }
    }

    break;
  }

  if (!microcanonicalParameterGridStructure_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "microcanonicalParameterGridStructure",
      "");
  }

  if (!singleChargeConfigurationArchiveList_.present ())
  {
    throw ::xsd::cxx::tree::expected_element< char > (
      "singleChargeConfigurationArchiveList",
      "");
  }
}

PartitionFunctionArchiveFile* PartitionFunctionArchiveFile::
_clone (::xml_schema::flags f,
        ::xml_schema::container* c) const
{
  return new class PartitionFunctionArchiveFile (*this, f, c);
}

PartitionFunctionArchiveFile& PartitionFunctionArchiveFile::
operator= (const PartitionFunctionArchiveFile& x)
{
  if (this != &x)
  {
    static_cast< ::xml_schema::type& > (*this) = x;
    this->microcanonicalParameterGridStructure_ = x.microcanonicalParameterGridStructure_;
    this->singleChargeConfigurationArchiveList_ = x.singleChargeConfigurationArchiveList_;
  }

  return *this;
}

PartitionFunctionArchiveFile::
~PartitionFunctionArchiveFile ()
{
}

#include <istream>
#include <xsd/cxx/xml/sax/std-input-source.hxx>
#include <xsd/cxx/tree/error-handler.hxx>

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (const ::std::string& u,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::tree::error_handler< char > h;

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      u, h, p, f));

  h.throw_if_failed< ::xsd::cxx::tree::parsing< char > > ();

  return ::std::unique_ptr< ::PartitionFunctionArchiveFile > (
    ::PartitionFunctionArchiveFile_ (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (const ::std::string& u,
                               ::xml_schema::error_handler& h,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      u, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::unique_ptr< ::PartitionFunctionArchiveFile > (
    ::PartitionFunctionArchiveFile_ (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (const ::std::string& u,
                               ::xercesc::DOMErrorHandler& h,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      u, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::unique_ptr< ::PartitionFunctionArchiveFile > (
    ::PartitionFunctionArchiveFile_ (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (::std::istream& is,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is);
  return ::PartitionFunctionArchiveFile_ (isrc, f, p);
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (::std::istream& is,
                               ::xml_schema::error_handler& h,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is);
  return ::PartitionFunctionArchiveFile_ (isrc, h, f, p);
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (::std::istream& is,
                               ::xercesc::DOMErrorHandler& h,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::sax::std_input_source isrc (is);
  return ::PartitionFunctionArchiveFile_ (isrc, h, f, p);
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (::std::istream& is,
                               const ::std::string& sid,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
  return ::PartitionFunctionArchiveFile_ (isrc, f, p);
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (::std::istream& is,
                               const ::std::string& sid,
                               ::xml_schema::error_handler& h,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0,
    (f & ::xml_schema::flags::keep_dom) == 0);

  ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
  return ::PartitionFunctionArchiveFile_ (isrc, h, f, p);
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (::std::istream& is,
                               const ::std::string& sid,
                               ::xercesc::DOMErrorHandler& h,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  ::xsd::cxx::xml::sax::std_input_source isrc (is, sid);
  return ::PartitionFunctionArchiveFile_ (isrc, h, f, p);
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (::xercesc::InputSource& i,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  ::xsd::cxx::tree::error_handler< char > h;

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      i, h, p, f));

  h.throw_if_failed< ::xsd::cxx::tree::parsing< char > > ();

  return ::std::unique_ptr< ::PartitionFunctionArchiveFile > (
    ::PartitionFunctionArchiveFile_ (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (::xercesc::InputSource& i,
                               ::xml_schema::error_handler& h,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      i, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::unique_ptr< ::PartitionFunctionArchiveFile > (
    ::PartitionFunctionArchiveFile_ (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (::xercesc::InputSource& i,
                               ::xercesc::DOMErrorHandler& h,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::parse< char > (
      i, h, p, f));

  if (!d.get ())
    throw ::xsd::cxx::tree::parsing< char > ();

  return ::std::unique_ptr< ::PartitionFunctionArchiveFile > (
    ::PartitionFunctionArchiveFile_ (
      std::move (d), f | ::xml_schema::flags::own_dom, p));
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (const ::xercesc::DOMDocument& doc,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties& p)
{
  if (f & ::xml_schema::flags::keep_dom)
  {
    ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
      static_cast< ::xercesc::DOMDocument* > (doc.cloneNode (true)));

    return ::std::unique_ptr< ::PartitionFunctionArchiveFile > (
      ::PartitionFunctionArchiveFile_ (
        std::move (d), f | ::xml_schema::flags::own_dom, p));
  }

  const ::xercesc::DOMElement& e (*doc.getDocumentElement ());
  const ::xsd::cxx::xml::qualified_name< char > n (
    ::xsd::cxx::xml::dom::name< char > (e));

  if (n.name () == "PartitionFunctionArchiveFile" &&
      n.namespace_ () == "")
  {
    ::std::unique_ptr< ::PartitionFunctionArchiveFile > r (
      ::xsd::cxx::tree::traits< ::PartitionFunctionArchiveFile, char >::create (
        e, f, 0));
    return r;
  }

  throw ::xsd::cxx::tree::unexpected_element < char > (
    n.name (),
    n.namespace_ (),
    "PartitionFunctionArchiveFile",
    "");
}

::std::unique_ptr< ::PartitionFunctionArchiveFile >
PartitionFunctionArchiveFile_ (::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d,
                               ::xml_schema::flags f,
                               const ::xml_schema::properties&)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > c (
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

  if (n.name () == "PartitionFunctionArchiveFile" &&
      n.namespace_ () == "")
  {
    ::std::unique_ptr< ::PartitionFunctionArchiveFile > r (
      ::xsd::cxx::tree::traits< ::PartitionFunctionArchiveFile, char >::create (
        e, f, 0));
    return r;
  }

  throw ::xsd::cxx::tree::unexpected_element < char > (
    n.name (),
    n.namespace_ (),
    "PartitionFunctionArchiveFile",
    "");
}

#include <ostream>
#include <xsd/cxx/tree/error-handler.hxx>
#include <xsd/cxx/xml/dom/serialization-source.hxx>

void
operator<< (::xercesc::DOMElement& e, const SingleChargeConfigurationArchive& i)
{
  e << static_cast< const ::xml_schema::type& > (i);

  // chargeConfiguration
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "chargeConfiguration",
        e));

    s << i.chargeConfiguration ();
  }

  // folderName
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "folderName",
        e));

    s << i.folderName ();
  }
}

void
operator<< (::xercesc::DOMElement& e, const SingleChargeConfigurationArchiveList& i)
{
  e << static_cast< const ::xml_schema::type& > (i);

  // singleChargeConfigurationArchive
  //
  for (SingleChargeConfigurationArchiveList::singleChargeConfigurationArchive_const_iterator
       b (i.singleChargeConfigurationArchive ().begin ()), n (i.singleChargeConfigurationArchive ().end ());
       b != n; ++b)
  {
    const SingleChargeConfigurationArchiveList::singleChargeConfigurationArchive_type& x (*b);

    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "singleChargeConfigurationArchive",
        e));

    s << x;
  }

  // length
  //
  {
    ::xercesc::DOMAttr& a (
      ::xsd::cxx::xml::dom::create_attribute (
        "length",
        e));

    a << i.length ();
  }
}

void
PartitionFunctionArchiveFile_ (::std::ostream& o,
                               const ::PartitionFunctionArchiveFile& s,
                               const ::xml_schema::namespace_infomap& m,
                               const ::std::string& e,
                               ::xml_schema::flags f)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0);

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::PartitionFunctionArchiveFile_ (s, m, f));

  ::xsd::cxx::tree::error_handler< char > h;

  ::xsd::cxx::xml::dom::ostream_format_target t (o);
  if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))
  {
    h.throw_if_failed< ::xsd::cxx::tree::serialization< char > > ();
  }
}

void
PartitionFunctionArchiveFile_ (::std::ostream& o,
                               const ::PartitionFunctionArchiveFile& s,
                               ::xml_schema::error_handler& h,
                               const ::xml_schema::namespace_infomap& m,
                               const ::std::string& e,
                               ::xml_schema::flags f)
{
  ::xsd::cxx::xml::auto_initializer i (
    (f & ::xml_schema::flags::dont_initialize) == 0);

  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::PartitionFunctionArchiveFile_ (s, m, f));
  ::xsd::cxx::xml::dom::ostream_format_target t (o);
  if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))
  {
    throw ::xsd::cxx::tree::serialization< char > ();
  }
}

void
PartitionFunctionArchiveFile_ (::std::ostream& o,
                               const ::PartitionFunctionArchiveFile& s,
                               ::xercesc::DOMErrorHandler& h,
                               const ::xml_schema::namespace_infomap& m,
                               const ::std::string& e,
                               ::xml_schema::flags f)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::PartitionFunctionArchiveFile_ (s, m, f));
  ::xsd::cxx::xml::dom::ostream_format_target t (o);
  if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))
  {
    throw ::xsd::cxx::tree::serialization< char > ();
  }
}

void
PartitionFunctionArchiveFile_ (::xercesc::XMLFormatTarget& t,
                               const ::PartitionFunctionArchiveFile& s,
                               const ::xml_schema::namespace_infomap& m,
                               const ::std::string& e,
                               ::xml_schema::flags f)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::PartitionFunctionArchiveFile_ (s, m, f));

  ::xsd::cxx::tree::error_handler< char > h;

  if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))
  {
    h.throw_if_failed< ::xsd::cxx::tree::serialization< char > > ();
  }
}

void
PartitionFunctionArchiveFile_ (::xercesc::XMLFormatTarget& t,
                               const ::PartitionFunctionArchiveFile& s,
                               ::xml_schema::error_handler& h,
                               const ::xml_schema::namespace_infomap& m,
                               const ::std::string& e,
                               ::xml_schema::flags f)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::PartitionFunctionArchiveFile_ (s, m, f));
  if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))
  {
    throw ::xsd::cxx::tree::serialization< char > ();
  }
}

void
PartitionFunctionArchiveFile_ (::xercesc::XMLFormatTarget& t,
                               const ::PartitionFunctionArchiveFile& s,
                               ::xercesc::DOMErrorHandler& h,
                               const ::xml_schema::namespace_infomap& m,
                               const ::std::string& e,
                               ::xml_schema::flags f)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::PartitionFunctionArchiveFile_ (s, m, f));
  if (!::xsd::cxx::xml::dom::serialize (t, *d, e, h, f))
  {
    throw ::xsd::cxx::tree::serialization< char > ();
  }
}

void
PartitionFunctionArchiveFile_ (::xercesc::DOMDocument& d,
                               const ::PartitionFunctionArchiveFile& s,
                               ::xml_schema::flags)
{
  ::xercesc::DOMElement& e (*d.getDocumentElement ());
  const ::xsd::cxx::xml::qualified_name< char > n (
    ::xsd::cxx::xml::dom::name< char > (e));

  if (n.name () == "PartitionFunctionArchiveFile" &&
      n.namespace_ () == "")
  {
    e << s;
  }
  else
  {
    throw ::xsd::cxx::tree::unexpected_element < char > (
      n.name (),
      n.namespace_ (),
      "PartitionFunctionArchiveFile",
      "");
  }
}

::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument >
PartitionFunctionArchiveFile_ (const ::PartitionFunctionArchiveFile& s,
                               const ::xml_schema::namespace_infomap& m,
                               ::xml_schema::flags f)
{
  ::xml_schema::dom::unique_ptr< ::xercesc::DOMDocument > d (
    ::xsd::cxx::xml::dom::serialize< char > (
      "PartitionFunctionArchiveFile",
      "",
      m, f));

  ::PartitionFunctionArchiveFile_ (*d, s, f);
  return d;
}

void
operator<< (::xercesc::DOMElement& e, const PartitionFunctionArchiveFile& i)
{
  e << static_cast< const ::xml_schema::type& > (i);

  // microcanonicalParameterGridStructure
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "microcanonicalParameterGridStructure",
        e));

    s << i.microcanonicalParameterGridStructure ();
  }

  // singleChargeConfigurationArchiveList
  //
  {
    ::xercesc::DOMElement& s (
      ::xsd::cxx::xml::dom::create_element (
        "singleChargeConfigurationArchiveList",
        e));

    s << i.singleChargeConfigurationArchiveList ();
  }
}

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

