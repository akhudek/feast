#define BOOST_TEST_MODULE extension_store test
#include <boost/test/included/unit_test.hpp>
#include <iostream>

#include "extension_store.h"

BOOST_AUTO_TEST_CASE( all_tests ) {
	using namespace std;
	extension_store store;

	extension_store::extension_ptr extA( new extension_store::extension() );
	extA->range.a = 10;
	extA->range.b = 12;
	extA->add_row( closed_interval(5,10), 6 );
	extA->add_row( closed_interval(7,12), 8 );
	extA->add_row( closed_interval(7,10), 6 );
	store.add( extA );

	// Test basic overlap code.
	try {
		store.get_bounds( point(11,11), true );
		BOOST_FAIL( "Overlap check failed." );
	} catch( point_covered_exception &e ) {}

	extension_store::bounds b = store.get_bounds( point(20,11), true );
	closed_interval test_int(8,20);
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 13 && test_int.b == 20 ) ;

	b = store.get_bounds( point(5,11), true );
	test_int.a = 0; test_int.b = 10;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.b == 6 && test_int.a == 0 );


	// Add new extension that overlaps extA completely.
	extension_store::extension_ptr extB( new extension_store::extension() );
	extB->range.a = 10;
	extB->range.b = 12;
	extB->add_row( closed_interval(15,20), 16 );
	extB->add_row( closed_interval(18,25), 19 );
	extB->add_row( closed_interval(18,22), 29 );
	store.add( extB );

	b = store.get_bounds( point(14,11), true );
	test_int.a = 8; test_int.b = 20;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 13 && test_int.b == 17 );

	// Add an extension that overlaps extA on the top exactly.
	extension_store::extension_ptr extC( new extension_store::extension() );
	extC->range.a = 10;
	extC->range.b = 15;
	extC->add_row( closed_interval(30,40), 36 );
	extC->add_row( closed_interval(32,45), 39 );
	extC->add_row( closed_interval(33,46), 40 );
	extC->add_row( closed_interval(33,47), 40 );
	extC->add_row( closed_interval(35,50), 43 );
	extC->add_row( closed_interval(36,45), 44 );
	store.add( extC );

	b = store.get_bounds( point(27,11), true );
	test_int.a = 20; test_int.b = 35;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 26 && test_int.b == 31 );

	b = store.get_bounds( point(27,14), true );
	test_int.a = 20; test_int.b = 40;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 20 && test_int.b == 34);

	// Break top extension and add a further extension on the bottom.
	extension_store::extension_ptr extD( new extension_store::extension() );
	extD->range.a = 12;
	extD->range.b = 20;
	extD->add_row( closed_interval(1,2), 1 );
	extD->add_row( closed_interval(1,3), 2 );
	extD->add_row( closed_interval(1,3), 2 ); //
	extD->add_row( closed_interval(1,3), 2 );
	extD->add_row( closed_interval(1,4), 2 );
	extD->add_row( closed_interval(1,4), 2 ); //
	extD->add_row( closed_interval(1,4), 2 );
	extD->add_row( closed_interval(1,4), 2 );
	extD->add_row( closed_interval(1,4), 2 ); //
	store.add( extD );

	b = store.get_bounds( point(5,11), true );
	test_int.a = 0; test_int.b = 40;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 0 && test_int.b == 6 );

	b = store.get_bounds( point(5,12), true );
	test_int.a = 0; test_int.b = 40;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 3 && test_int.b == 6 );

	b = store.get_bounds( point(5,14), true );
	test_int.a = 0; test_int.b = 40;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 4 && test_int.b == 34 );

	b = store.get_bounds( point(5,18), true );
	test_int.a = 0; test_int.b = 40;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 5 && test_int.b == 40 );

	// Add overlap where top has overhang.
	extension_store::extension_ptr extE( new extension_store::extension() );
	extE->range.a = 6;
	extE->range.b = 11;
	extE->add_row( closed_interval(80,90), 85 );
	extE->add_row( closed_interval(80,90), 85 );
	extE->add_row( closed_interval(80,90), 85 );
	extE->add_row( closed_interval(80,90), 85 );//
	extE->add_row( closed_interval(82,92), 85 );
	extE->add_row( closed_interval(82,92), 85 );
	store.add( extE );

	b = store.get_bounds( point(60,8), true );
	test_int.a = 0; test_int.b = 100;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 0 && test_int.b == 79 );

	b = store.get_bounds( point(60,10), true );
	test_int.a = 0; test_int.b = 100;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 41 && test_int.b == 81 );

	// Need to test iteration.
	closed_interval bint(25,25);
	b = store.get_bounds( point(25,0), true );
	test_int.a = 0; test_int.b = 100;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 0 && test_int.b == 100 );

	for( int i = 0; i < 6; ++i ) b.step(bint);
	test_int.a = 0; test_int.b = 100;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 0 && test_int.b == 79 );

	for( int i = 0; i < 4; ++i ) b.step(bint);
	test_int.a = 0; test_int.b = 100;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 21 && test_int.b == 29 );

	for( int i = 0; i < 2; ++i ) b.step(bint);
	test_int.a = 0; test_int.b = 100;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 23 && test_int.b == 32 );

	for( int i = 0; i < 100; ++i ) b.step(bint);
	test_int.a = 0; test_int.b = 100;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 0 && test_int.b == 100 );

	// Test reverse step.
	b = store.get_bounds( point(25,25), false );
	test_int.a = 0; test_int.b = 100;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 0 && test_int.b == 100 );

	for( int i = 0; i < 5; ++i ) b.step(bint);
	test_int.a = 0; test_int.b = 100;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 5 && test_int.b == 100 );

	for( int i = 0; i < 6; ++i ) b.step(bint);
	test_int.a = 0; test_int.b = 100;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 4 && test_int.b == 34 );

	for( int i = 0; i < 13; ++i ) b.step(bint);
	test_int.a = 0; test_int.b = 100;
	test_int = b.clip(test_int);
	BOOST_CHECK( test_int.a == 0 && test_int.b == 100 );

	store.clear();
	extension_store::extension_ptr extF( new extension_store::extension() );
	extF->range.a = 10;
	extF->range.b = 20;
	extF->add_row( closed_interval(5,10), 6 );
	extF->add_row( closed_interval(5,10), 6 );
	extF->add_row( closed_interval(5,10), 6 );
	extF->add_row( closed_interval(5,10), 6 );
	extF->add_row( closed_interval(5,10), 6 );
	extF->add_row( closed_interval(5,10), 6 );
	extF->add_row( closed_interval(5,10), 6 );
	extF->add_row( closed_interval(5,10), 6 );
	extF->add_row( closed_interval(5,10), 6 );
	extF->add_row( closed_interval(5,10), 6 );
	extF->add_row( closed_interval(5,10), 6 );
	store.add( extF );

	extension_store::extension_ptr extH( new extension_store::extension() );
	extH->range.a = 10;
	extH->range.b = 20;
	extH->add_row( closed_interval(11,13), 11 );
	extH->add_row( closed_interval(11,13), 11 );
	extH->add_row( closed_interval(11,13), 11 );
	extH->add_row( closed_interval(11,13), 11 );
	extH->add_row( closed_interval(11,13), 11 );
	extH->add_row( closed_interval(11,13), 11 );
	extH->add_row( closed_interval(11,13), 11 );
	extH->add_row( closed_interval(11,13), 11 );
	extH->add_row( closed_interval(11,13), 11 );
	extH->add_row( closed_interval(11,13), 11 );
	extH->add_row( closed_interval(11,13), 11 );
	store.add( extH );

	extension_store::extension_ptr extG( new extension_store::extension() );
	extG->range.a = 13;
	extG->range.b = 16;
	extG->add_row( closed_interval(15,20), 17 );
	extG->add_row( closed_interval(15,20), 17 );
	extG->add_row( closed_interval(15,20), 17 );
	extG->add_row( closed_interval(15,20), 17 );
	store.add( extG );
	store.check();

}
