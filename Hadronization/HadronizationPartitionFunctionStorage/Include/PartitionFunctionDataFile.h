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

#ifndef PARTITION_FUNCTION_DATA_FILE_H
#define PARTITION_FUNCTION_DATA_FILE_H

// Begin prologue.
//
//
// End prologue.

#include <xsd/cxx/config.hxx>

#if (XSD_INT_VERSION != 3030000L)
#error XSD runtime version mismatch
#endif

#include <xsd/cxx/pre.hxx>

#ifndef XSD_USE_CHAR
#define XSD_USE_CHAR
#endif

#ifndef XSD_CXX_TREE_USE_CHAR
#define XSD_CXX_TREE_USE_CHAR
#endif

#include <xsd/cxx/xml/char-utf8.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/types.hxx>

#include <xsd/cxx/xml/error-handler.hxx>

#include <xsd/cxx/xml/dom/auto-ptr.hxx>

#include <xsd/cxx/tree/parsing.hxx>
#include <xsd/cxx/tree/parsing/byte.hxx>
#include <xsd/cxx/tree/parsing/unsigned-byte.hxx>
#include <xsd/cxx/tree/parsing/short.hxx>
#include <xsd/cxx/tree/parsing/unsigned-short.hxx>
#include <xsd/cxx/tree/parsing/int.hxx>
#include <xsd/cxx/tree/parsing/unsigned-int.hxx>
#include <xsd/cxx/tree/parsing/long.hxx>
#include <xsd/cxx/tree/parsing/unsigned-long.hxx>
#include <xsd/cxx/tree/parsing/boolean.hxx>
#include <xsd/cxx/tree/parsing/float.hxx>
#include <xsd/cxx/tree/parsing/double.hxx>
#include <xsd/cxx/tree/parsing/decimal.hxx>

#include <xsd/cxx/xml/dom/serialization-header.hxx>
#include <xsd/cxx/tree/serialization.hxx>
#include <xsd/cxx/tree/serialization/byte.hxx>
#include <xsd/cxx/tree/serialization/unsigned-byte.hxx>
#include <xsd/cxx/tree/serialization/short.hxx>
#include <xsd/cxx/tree/serialization/unsigned-short.hxx>
#include <xsd/cxx/tree/serialization/int.hxx>
#include <xsd/cxx/tree/serialization/unsigned-int.hxx>
#include <xsd/cxx/tree/serialization/long.hxx>
#include <xsd/cxx/tree/serialization/unsigned-long.hxx>
#include <xsd/cxx/tree/serialization/boolean.hxx>
#include <xsd/cxx/tree/serialization/float.hxx>
#include <xsd/cxx/tree/serialization/double.hxx>
#include <xsd/cxx/tree/serialization/decimal.hxx>

namespace xml_schema
{
  // anyType and anySimpleType.
  //
  typedef ::xsd::cxx::tree::type type;
  typedef ::xsd::cxx::tree::simple_type< type > simple_type;
  typedef ::xsd::cxx::tree::type container;

  // 8-bit
  //
  typedef signed char byte;
  typedef unsigned char unsigned_byte;

  // 16-bit
  //
  typedef short short_;
  typedef unsigned short unsigned_short;

  // 32-bit
  //
  typedef int int_;
  typedef unsigned int unsigned_int;

  // 64-bit
  //
  typedef long long long_;
  typedef unsigned long long unsigned_long;

  // Supposed to be arbitrary-length integral types.
  //
  typedef long long integer;
  typedef long long non_positive_integer;
  typedef unsigned long long non_negative_integer;
  typedef unsigned long long positive_integer;
  typedef long long negative_integer;

  // Boolean.
  //
  typedef bool boolean;

  // Floating-point types.
  //
  typedef float float_;
  typedef double double_;
  typedef double decimal;

  // String types.
  //
  typedef ::xsd::cxx::tree::string< char, simple_type > string;
  typedef ::xsd::cxx::tree::normalized_string< char, string > normalized_string;
  typedef ::xsd::cxx::tree::token< char, normalized_string > token;
  typedef ::xsd::cxx::tree::name< char, token > name;
  typedef ::xsd::cxx::tree::nmtoken< char, token > nmtoken;
  typedef ::xsd::cxx::tree::nmtokens< char, simple_type, nmtoken > nmtokens;
  typedef ::xsd::cxx::tree::ncname< char, name > ncname;
  typedef ::xsd::cxx::tree::language< char, token > language;

  // ID/IDREF.
  //
  typedef ::xsd::cxx::tree::id< char, ncname > id;
  typedef ::xsd::cxx::tree::idref< char, ncname, type > idref;
  typedef ::xsd::cxx::tree::idrefs< char, simple_type, idref > idrefs;

  // URI.
  //
  typedef ::xsd::cxx::tree::uri< char, simple_type > uri;

  // Qualified name.
  //
  typedef ::xsd::cxx::tree::qname< char, simple_type, uri, ncname > qname;

  // Binary.
  //
  typedef ::xsd::cxx::tree::buffer< char > buffer;
  typedef ::xsd::cxx::tree::base64_binary< char, simple_type > base64_binary;
  typedef ::xsd::cxx::tree::hex_binary< char, simple_type > hex_binary;

  // Date/time.
  //
  typedef ::xsd::cxx::tree::time_zone time_zone;
  typedef ::xsd::cxx::tree::date< char, simple_type > date;
  typedef ::xsd::cxx::tree::date_time< char, simple_type > date_time;
  typedef ::xsd::cxx::tree::duration< char, simple_type > duration;
  typedef ::xsd::cxx::tree::gday< char, simple_type > gday;
  typedef ::xsd::cxx::tree::gmonth< char, simple_type > gmonth;
  typedef ::xsd::cxx::tree::gmonth_day< char, simple_type > gmonth_day;
  typedef ::xsd::cxx::tree::gyear< char, simple_type > gyear;
  typedef ::xsd::cxx::tree::gyear_month< char, simple_type > gyear_month;
  typedef ::xsd::cxx::tree::time< char, simple_type > time;

  // Entity.
  //
  typedef ::xsd::cxx::tree::entity< char, ncname > entity;
  typedef ::xsd::cxx::tree::entities< char, simple_type, entity > entities;

  // Namespace information and list stream. Used in
  // serialization functions.
  //
  typedef ::xsd::cxx::xml::dom::namespace_info< char > namespace_info;
  typedef ::xsd::cxx::xml::dom::namespace_infomap< char > namespace_infomap;
  typedef ::xsd::cxx::tree::list_stream< char > list_stream;
  typedef ::xsd::cxx::tree::as_double< double_ > as_double;
  typedef ::xsd::cxx::tree::as_decimal< decimal > as_decimal;
  typedef ::xsd::cxx::tree::facet facet;

  // Flags and properties.
  //
  typedef ::xsd::cxx::tree::flags flags;
  typedef ::xsd::cxx::tree::properties< char > properties;

  // Parsing/serialization diagnostics.
  //
  typedef ::xsd::cxx::tree::severity severity;
  typedef ::xsd::cxx::tree::error< char > error;
  typedef ::xsd::cxx::tree::diagnostics< char > diagnostics;

  // Exceptions.
  //
  typedef ::xsd::cxx::tree::exception< char > exception;
  typedef ::xsd::cxx::tree::bounds< char > bounds;
  typedef ::xsd::cxx::tree::duplicate_id< char > duplicate_id;
  typedef ::xsd::cxx::tree::parsing< char > parsing;
  typedef ::xsd::cxx::tree::expected_element< char > expected_element;
  typedef ::xsd::cxx::tree::unexpected_element< char > unexpected_element;
  typedef ::xsd::cxx::tree::expected_attribute< char > expected_attribute;
  typedef ::xsd::cxx::tree::unexpected_enumerator< char > unexpected_enumerator;
  typedef ::xsd::cxx::tree::expected_text_content< char > expected_text_content;
  typedef ::xsd::cxx::tree::no_prefix_mapping< char > no_prefix_mapping;
  typedef ::xsd::cxx::tree::serialization< char > serialization;

  // Error handler callback interface.
  //
  typedef ::xsd::cxx::xml::error_handler< char > error_handler;

  // DOM interaction.
  //
  namespace dom
  {
    // Automatic pointer for DOMDocument.
    //
    using ::xsd::cxx::xml::dom::auto_ptr;

#ifndef XSD_CXX_TREE_TREE_NODE_KEY__XML_SCHEMA
#define XSD_CXX_TREE_TREE_NODE_KEY__XML_SCHEMA
    // DOM user data key for back pointers to tree nodes.
    //
    const XMLCh* const tree_node_key = ::xsd::cxx::tree::user_data_keys::node;
#endif
  }
}

// Forward declarations.
//
class ImportanceSamplingFunctionParameters;
class PartitionFunctionGridElement;
class PartitionFunctionDataList;
class PartitionFunctionDataFile;

#include <memory>    // std::auto_ptr
#include <limits>    // std::numeric_limits
#include <algorithm> // std::binary_search

#include <xsd/cxx/xml/char-utf8.hxx>

#include <xsd/cxx/tree/exceptions.hxx>
#include <xsd/cxx/tree/elements.hxx>
#include <xsd/cxx/tree/containers.hxx>
#include <xsd/cxx/tree/list.hxx>

#include <xsd/cxx/xml/dom/parsing-header.hxx>

class ImportanceSamplingFunctionParameters: public ::xml_schema::type
{
  public:
  // electricChargeFugacity
  // 
  typedef ::xml_schema::double_ electricChargeFugacity_type;
  typedef ::xsd::cxx::tree::traits< electricChargeFugacity_type, char, ::xsd::cxx::tree::schema_type::double_ > electricChargeFugacity_traits;

  const electricChargeFugacity_type&
  electricChargeFugacity () const;

  electricChargeFugacity_type&
  electricChargeFugacity ();

  void
  electricChargeFugacity (const electricChargeFugacity_type& x);

  // baryonicChargeFugacity
  // 
  typedef ::xml_schema::double_ baryonicChargeFugacity_type;
  typedef ::xsd::cxx::tree::traits< baryonicChargeFugacity_type, char, ::xsd::cxx::tree::schema_type::double_ > baryonicChargeFugacity_traits;

  const baryonicChargeFugacity_type&
  baryonicChargeFugacity () const;

  baryonicChargeFugacity_type&
  baryonicChargeFugacity ();

  void
  baryonicChargeFugacity (const baryonicChargeFugacity_type& x);

  // strangeChargeFugacity
  // 
  typedef ::xml_schema::double_ strangeChargeFugacity_type;
  typedef ::xsd::cxx::tree::traits< strangeChargeFugacity_type, char, ::xsd::cxx::tree::schema_type::double_ > strangeChargeFugacity_traits;

  const strangeChargeFugacity_type&
  strangeChargeFugacity () const;

  strangeChargeFugacity_type&
  strangeChargeFugacity ();

  void
  strangeChargeFugacity (const strangeChargeFugacity_type& x);

  // charmChargeFugacity
  // 
  typedef ::xml_schema::double_ charmChargeFugacity_type;
  typedef ::xsd::cxx::tree::traits< charmChargeFugacity_type, char, ::xsd::cxx::tree::schema_type::double_ > charmChargeFugacity_traits;

  const charmChargeFugacity_type&
  charmChargeFugacity () const;

  charmChargeFugacity_type&
  charmChargeFugacity ();

  void
  charmChargeFugacity (const charmChargeFugacity_type& x);

  // bottomChargeFugacity
  // 
  typedef ::xml_schema::double_ bottomChargeFugacity_type;
  typedef ::xsd::cxx::tree::traits< bottomChargeFugacity_type, char, ::xsd::cxx::tree::schema_type::double_ > bottomChargeFugacity_traits;

  const bottomChargeFugacity_type&
  bottomChargeFugacity () const;

  bottomChargeFugacity_type&
  bottomChargeFugacity ();

  void
  bottomChargeFugacity (const bottomChargeFugacity_type& x);

  // temperature
  // 
  typedef ::xml_schema::double_ temperature_type;
  typedef ::xsd::cxx::tree::traits< temperature_type, char, ::xsd::cxx::tree::schema_type::double_ > temperature_traits;

  const temperature_type&
  temperature () const;

  temperature_type&
  temperature ();

  void
  temperature (const temperature_type& x);

  // Constructors.
  //
  ImportanceSamplingFunctionParameters (const electricChargeFugacity_type&,
                                        const baryonicChargeFugacity_type&,
                                        const strangeChargeFugacity_type&,
                                        const charmChargeFugacity_type&,
                                        const bottomChargeFugacity_type&,
                                        const temperature_type&);

  ImportanceSamplingFunctionParameters (const ::xercesc::DOMElement& e,
                                        ::xml_schema::flags f = 0,
                                        ::xml_schema::container* c = 0);

  ImportanceSamplingFunctionParameters (const ImportanceSamplingFunctionParameters& x,
                                        ::xml_schema::flags f = 0,
                                        ::xml_schema::container* c = 0);

  virtual ImportanceSamplingFunctionParameters*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  virtual 
  ~ImportanceSamplingFunctionParameters ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< electricChargeFugacity_type > electricChargeFugacity_;
  ::xsd::cxx::tree::one< baryonicChargeFugacity_type > baryonicChargeFugacity_;
  ::xsd::cxx::tree::one< strangeChargeFugacity_type > strangeChargeFugacity_;
  ::xsd::cxx::tree::one< charmChargeFugacity_type > charmChargeFugacity_;
  ::xsd::cxx::tree::one< bottomChargeFugacity_type > bottomChargeFugacity_;
  ::xsd::cxx::tree::one< temperature_type > temperature_;
};

class PartitionFunctionGridElement: public ::xml_schema::type
{
  public:
  // partitionFunctionGridElementIndex
  // 
  typedef ::xml_schema::unsigned_int partitionFunctionGridElementIndex_type;
  typedef ::xsd::cxx::tree::traits< partitionFunctionGridElementIndex_type, char > partitionFunctionGridElementIndex_traits;

  const partitionFunctionGridElementIndex_type&
  partitionFunctionGridElementIndex () const;

  partitionFunctionGridElementIndex_type&
  partitionFunctionGridElementIndex ();

  void
  partitionFunctionGridElementIndex (const partitionFunctionGridElementIndex_type& x);

  // partitionFunctionValue
  // 
  typedef ::xml_schema::double_ partitionFunctionValue_type;
  typedef ::xsd::cxx::tree::traits< partitionFunctionValue_type, char, ::xsd::cxx::tree::schema_type::double_ > partitionFunctionValue_traits;

  const partitionFunctionValue_type&
  partitionFunctionValue () const;

  partitionFunctionValue_type&
  partitionFunctionValue ();

  void
  partitionFunctionValue (const partitionFunctionValue_type& x);

  // partitionFunctionErrorValue
  // 
  typedef ::xml_schema::double_ partitionFunctionErrorValue_type;
  typedef ::xsd::cxx::tree::traits< partitionFunctionErrorValue_type, char, ::xsd::cxx::tree::schema_type::double_ > partitionFunctionErrorValue_traits;

  const partitionFunctionErrorValue_type&
  partitionFunctionErrorValue () const;

  partitionFunctionErrorValue_type&
  partitionFunctionErrorValue ();

  void
  partitionFunctionErrorValue (const partitionFunctionErrorValue_type& x);

  // numberChannelSamplings
  // 
  typedef ::xml_schema::unsigned_long numberChannelSamplings_type;
  typedef ::xsd::cxx::tree::traits< numberChannelSamplings_type, char > numberChannelSamplings_traits;

  const numberChannelSamplings_type&
  numberChannelSamplings () const;

  numberChannelSamplings_type&
  numberChannelSamplings ();

  void
  numberChannelSamplings (const numberChannelSamplings_type& x);

  // importanceSamplingFunctionParameters
  // 
  typedef ::ImportanceSamplingFunctionParameters importanceSamplingFunctionParameters_type;
  typedef ::xsd::cxx::tree::optional< importanceSamplingFunctionParameters_type > importanceSamplingFunctionParameters_optional;
  typedef ::xsd::cxx::tree::traits< importanceSamplingFunctionParameters_type, char > importanceSamplingFunctionParameters_traits;

  const importanceSamplingFunctionParameters_optional&
  importanceSamplingFunctionParameters () const;

  importanceSamplingFunctionParameters_optional&
  importanceSamplingFunctionParameters ();

  void
  importanceSamplingFunctionParameters (const importanceSamplingFunctionParameters_type& x);

  void
  importanceSamplingFunctionParameters (const importanceSamplingFunctionParameters_optional& x);

  void
  importanceSamplingFunctionParameters (::std::auto_ptr< importanceSamplingFunctionParameters_type > p);

  // Constructors.
  //
  PartitionFunctionGridElement (const partitionFunctionGridElementIndex_type&,
                                const partitionFunctionValue_type&,
                                const partitionFunctionErrorValue_type&,
                                const numberChannelSamplings_type&);

  PartitionFunctionGridElement (const ::xercesc::DOMElement& e,
                                ::xml_schema::flags f = 0,
                                ::xml_schema::container* c = 0);

  PartitionFunctionGridElement (const PartitionFunctionGridElement& x,
                                ::xml_schema::flags f = 0,
                                ::xml_schema::container* c = 0);

  virtual PartitionFunctionGridElement*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  virtual 
  ~PartitionFunctionGridElement ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< partitionFunctionGridElementIndex_type > partitionFunctionGridElementIndex_;
  ::xsd::cxx::tree::one< partitionFunctionValue_type > partitionFunctionValue_;
  ::xsd::cxx::tree::one< partitionFunctionErrorValue_type > partitionFunctionErrorValue_;
  ::xsd::cxx::tree::one< numberChannelSamplings_type > numberChannelSamplings_;
  importanceSamplingFunctionParameters_optional importanceSamplingFunctionParameters_;
};

class PartitionFunctionDataList: public ::xml_schema::type
{
  public:
  // partitionFunctionData
  // 
  typedef ::PartitionFunctionGridElement partitionFunctionData_type;
  typedef ::xsd::cxx::tree::sequence< partitionFunctionData_type > partitionFunctionData_sequence;
  typedef partitionFunctionData_sequence::iterator partitionFunctionData_iterator;
  typedef partitionFunctionData_sequence::const_iterator partitionFunctionData_const_iterator;
  typedef ::xsd::cxx::tree::traits< partitionFunctionData_type, char > partitionFunctionData_traits;

  const partitionFunctionData_sequence&
  partitionFunctionData () const;

  partitionFunctionData_sequence&
  partitionFunctionData ();

  void
  partitionFunctionData (const partitionFunctionData_sequence& s);

  // length
  // 
  typedef ::xml_schema::unsigned_int length_type;
  typedef ::xsd::cxx::tree::traits< length_type, char > length_traits;

  const length_type&
  length () const;

  length_type&
  length ();

  void
  length (const length_type& x);

  // Constructors.
  //
  PartitionFunctionDataList (const length_type&);

  PartitionFunctionDataList (const ::xercesc::DOMElement& e,
                             ::xml_schema::flags f = 0,
                             ::xml_schema::container* c = 0);

  PartitionFunctionDataList (const PartitionFunctionDataList& x,
                             ::xml_schema::flags f = 0,
                             ::xml_schema::container* c = 0);

  virtual PartitionFunctionDataList*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  virtual 
  ~PartitionFunctionDataList ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  partitionFunctionData_sequence partitionFunctionData_;
  ::xsd::cxx::tree::one< length_type > length_;
};

class PartitionFunctionDataFile: public ::xml_schema::type
{
  public:
  // massValue
  // 
  typedef ::xml_schema::double_ massValue_type;
  typedef ::xsd::cxx::tree::traits< massValue_type, char, ::xsd::cxx::tree::schema_type::double_ > massValue_traits;

  const massValue_type&
  massValue () const;

  massValue_type&
  massValue ();

  void
  massValue (const massValue_type& x);

  // partitionFunctionDataList
  // 
  typedef ::PartitionFunctionDataList partitionFunctionDataList_type;
  typedef ::xsd::cxx::tree::traits< partitionFunctionDataList_type, char > partitionFunctionDataList_traits;

  const partitionFunctionDataList_type&
  partitionFunctionDataList () const;

  partitionFunctionDataList_type&
  partitionFunctionDataList ();

  void
  partitionFunctionDataList (const partitionFunctionDataList_type& x);

  void
  partitionFunctionDataList (::std::auto_ptr< partitionFunctionDataList_type > p);

  // Constructors.
  //
  PartitionFunctionDataFile (const massValue_type&,
                             const partitionFunctionDataList_type&);

  PartitionFunctionDataFile (const massValue_type&,
                             ::std::auto_ptr< partitionFunctionDataList_type >&);

  PartitionFunctionDataFile (const ::xercesc::DOMElement& e,
                             ::xml_schema::flags f = 0,
                             ::xml_schema::container* c = 0);

  PartitionFunctionDataFile (const PartitionFunctionDataFile& x,
                             ::xml_schema::flags f = 0,
                             ::xml_schema::container* c = 0);

  virtual PartitionFunctionDataFile*
  _clone (::xml_schema::flags f = 0,
          ::xml_schema::container* c = 0) const;

  virtual 
  ~PartitionFunctionDataFile ();

  // Implementation.
  //
  protected:
  void
  parse (::xsd::cxx::xml::dom::parser< char >&,
         ::xml_schema::flags);

  protected:
  ::xsd::cxx::tree::one< massValue_type > massValue_;
  ::xsd::cxx::tree::one< partitionFunctionDataList_type > partitionFunctionDataList_;
};

#include <iosfwd>

#include <xercesc/sax/InputSource.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>

// Parse a URI or a local file.
//

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (const ::std::string& uri,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (const ::std::string& uri,
                            ::xml_schema::error_handler& eh,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (const ::std::string& uri,
                            ::xercesc::DOMErrorHandler& eh,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse std::istream.
//

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (::std::istream& is,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (::std::istream& is,
                            ::xml_schema::error_handler& eh,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (::std::istream& is,
                            ::xercesc::DOMErrorHandler& eh,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (::std::istream& is,
                            const ::std::string& id,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (::std::istream& is,
                            const ::std::string& id,
                            ::xml_schema::error_handler& eh,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (::std::istream& is,
                            const ::std::string& id,
                            ::xercesc::DOMErrorHandler& eh,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse xercesc::InputSource.
//

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (::xercesc::InputSource& is,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (::xercesc::InputSource& is,
                            ::xml_schema::error_handler& eh,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (::xercesc::InputSource& is,
                            ::xercesc::DOMErrorHandler& eh,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

// Parse xercesc::DOMDocument.
//

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (const ::xercesc::DOMDocument& d,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

::std::auto_ptr< ::PartitionFunctionDataFile >
PartitionFunctionDataFile_ (::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument >& d,
                            ::xml_schema::flags f = 0,
                            const ::xml_schema::properties& p = ::xml_schema::properties ());

#include <iosfwd>

#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>
#include <xercesc/framework/XMLFormatter.hpp>

#include <xsd/cxx/xml/dom/auto-ptr.hxx>

void
operator<< (::xercesc::DOMElement&, const ImportanceSamplingFunctionParameters&);

void
operator<< (::xercesc::DOMElement&, const PartitionFunctionGridElement&);

void
operator<< (::xercesc::DOMElement&, const PartitionFunctionDataList&);

// Serialize to std::ostream.
//

void
PartitionFunctionDataFile_ (::std::ostream& os,
                            const ::PartitionFunctionDataFile& x, 
                            const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
                            const ::std::string& e = "UTF-8",
                            ::xml_schema::flags f = 0);

void
PartitionFunctionDataFile_ (::std::ostream& os,
                            const ::PartitionFunctionDataFile& x, 
                            ::xml_schema::error_handler& eh,
                            const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
                            const ::std::string& e = "UTF-8",
                            ::xml_schema::flags f = 0);

void
PartitionFunctionDataFile_ (::std::ostream& os,
                            const ::PartitionFunctionDataFile& x, 
                            ::xercesc::DOMErrorHandler& eh,
                            const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
                            const ::std::string& e = "UTF-8",
                            ::xml_schema::flags f = 0);

// Serialize to xercesc::XMLFormatTarget.
//

void
PartitionFunctionDataFile_ (::xercesc::XMLFormatTarget& ft,
                            const ::PartitionFunctionDataFile& x, 
                            const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
                            const ::std::string& e = "UTF-8",
                            ::xml_schema::flags f = 0);

void
PartitionFunctionDataFile_ (::xercesc::XMLFormatTarget& ft,
                            const ::PartitionFunctionDataFile& x, 
                            ::xml_schema::error_handler& eh,
                            const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
                            const ::std::string& e = "UTF-8",
                            ::xml_schema::flags f = 0);

void
PartitionFunctionDataFile_ (::xercesc::XMLFormatTarget& ft,
                            const ::PartitionFunctionDataFile& x, 
                            ::xercesc::DOMErrorHandler& eh,
                            const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
                            const ::std::string& e = "UTF-8",
                            ::xml_schema::flags f = 0);

// Serialize to an existing xercesc::DOMDocument.
//

void
PartitionFunctionDataFile_ (::xercesc::DOMDocument& d,
                            const ::PartitionFunctionDataFile& x,
                            ::xml_schema::flags f = 0);

// Serialize to a new xercesc::DOMDocument.
//

::xml_schema::dom::auto_ptr< ::xercesc::DOMDocument >
PartitionFunctionDataFile_ (const ::PartitionFunctionDataFile& x, 
                            const ::xml_schema::namespace_infomap& m = ::xml_schema::namespace_infomap (),
                            ::xml_schema::flags f = 0);

void
operator<< (::xercesc::DOMElement&, const PartitionFunctionDataFile&);

#include <xsd/cxx/post.hxx>

// Begin epilogue.
//
//
// End epilogue.

#endif // PARTITION_FUNCTION_DATA_FILE_H