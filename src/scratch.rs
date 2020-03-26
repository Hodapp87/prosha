use std::rc::Rc;

/*
struct R<'a> {
    b: &'a dyn Fn() -> R<'a>,
}

#[derive(Copy, Clone)]
struct Foo {}
impl<'a> Foo {
    // These are valid, but not especially useful (if I am
    // transferring ownership then I cannot have any branching):
    fn fn1(self) -> R<'a> {
        R { b: & move || self.fn1() }
    }
    fn fn2(self) -> R<'a> {
        R { b: &|| self.fn2() }
    }
}
*/

// Below (using box instead of a trait object) follows similar rules:
struct S<'a> {
    b: Box<dyn Fn() -> S<'a>>,
}
#[derive(Copy, Clone)]
struct Foo2 {}
impl<'a> Foo2 {
    fn fn1(self) -> S<'a> {
        S { b: Box::new(move || self.fn1()) }
    }
    // Not valid (error[E0373]: closure may outlive the current
    // function, but it borrows `self`, which is owned by the current
    // function):
    //fn fn2(self) -> S<'a> {
    //    S { b: Box::new(|| self.fn2()) }
    //}
    // Not valid:
    //fn fn3(&self) -> S<'a> {
    //    S { b: Box::new(move || self.fn3()) }
    //}
    // Not valid:
    //fn fn4(&self) -> S<'a> {
    //    S { b: Box::new(|| self.fn4()) }
    //}
}

struct T<'a> {
    b: Rc<dyn Fn() -> T<'a> + 'a>,
}
#[derive(Copy, Clone)]
struct Foo3 {}
impl<'a> Foo3 {
    fn fn1(self) -> T<'a> {
        T { b: Rc::new(move || self.fn1()) }
    }
    // Not valid (E0373):
    //fn fn2(self) -> T<'a> {
    //    T { b: Rc::new(|| self.fn2()) }
    //}
    // Not valid:
    //fn fn3(&self) -> T<'a> {
    //    T { b: Rc::new(move || self.fn3()) }
    //}
    // Not valid:
    //fn fn4(&self) -> T<'a> {
    //    T { b: Rc::new(|| self.fn4()) }
    //}
    // But this is now valid because T can be cloned:
    fn fn5(self) -> (T<'a>, T<'a>) {
        let p = Rc::new(move || self.fn1());
        let p2 = p.clone();
        (T { b: p }, T { b: p2 })
    }
}

// Further, this is now valid too (lifetimes removed):
struct U {
    b: Rc<dyn Fn() -> U>,
}
#[derive(Copy, Clone)]
struct Foo4 {}
impl Foo4 {
    fn fn1(self) -> U {
        U { b: Rc::new(move || self.fn1()) }
    }
    fn fn5(self) -> (U, U) {
        let p = Rc::new(move || self.fn1());
        let p2 = p.clone();
        (U { b: p }, U { b: p2 })
    }
}

// I can get rid of Copy/Clone if I use FnOnce:
struct V {
    b: Rc<dyn FnOnce() -> V>,
}
struct Foo5 {}
impl Foo5 {
    fn fn1(self) -> V {
        V { b: Rc::new(move || self.fn1()) }
    }
    fn fn2(self) -> (V, V) {
        let p = Rc::new(move || self.fn1());
        let p2 = p.clone();
        (V { b: p }, V { b: p2 })
    }
    // and then either kind is fine:
    fn fn3(self) -> V {
        V { b: Rc::new(|| self.fn3()) }
    }
    fn fn4(self) -> (V, V) {
        let p = Rc::new(|| self.fn3());
        let p2 = p.clone();
        (V { b: p }, V { b: p2 })
        // but this confuses me a bit. doesn't this then let me call
        // an FnOnce... more than once?
    }
}

// This is valid and I can recurse:
struct W {
    b: Box<dyn Fn() -> W>,
}
struct Foo6 {}
impl Foo6 {
    fn fn1(s: &Rc<Self>) -> W {
        let s2 = Rc::clone(&s);
        W { b: Box::new(move || Self::fn1(&s2)) }
    }
    fn fn2(s: &Rc<Self>) -> (W, W) {
        let s2 = Rc::clone(&s);
        let w2 = W { b: Box::new(move || Self::fn1(&s2)) };
        let s3 = Rc::clone(&s);
        let w3 = W { b: Box::new(move || Self::fn1(&s3)) };
        (w2, w3)
    }
}

fn foo6() {

    // Whatever (note that it doesn't automatically do Copy):
    struct State {
        v: u32,
    }

    // Purposely put state somewhere it goes out of scope:
    let s = {
        let s_orig = State {
            v: 105,
        };
        Rc::new(s_orig)
    };
    /*
    let fn1 = |f: &dyn Fn(&dyn Fn() -> W) -> (&dyn Fn() -> W)| -> (&dyn Fn() -> W) {
        &(|| -> W {
            let s2 = Rc::clone(&s);
            W { b: Box::new(move || f(f)) }
        })
    };

    let f2 = fn1(fn1);
    */
}

fn foo7(t: impl Clone) -> impl Clone {
    t.clone()
}

fn foo7b<T: Clone>(t: T) -> T {
    t.clone()
}

fn foo7c<T>(t: T) -> T where T: Clone {
    t.clone()
}

// A simple implementation of the Y Combinator
// λf.(λx.xx)(λx.f(xx))
// <=> λf.(λx.f(xx))(λx.f(xx))
 
// CREDITS: A better version of the previous code that was posted here, with detailed explanation.
// See <y> and also <y_apply>.
 
// A function type that takes its own type as an input is an infinite recursive type.
// We introduce a trait that will allow us to have an input with the same type as self, and break the recursion.
// The input is going to be a trait object that implements the desired function in the interface.
// NOTE: We will be coercing a reference to a closure into this trait object.
 
trait Apply<T, R> {
    fn apply(&self, f: &dyn Apply<T, R>, t: T) -> R;
}
 
// In Rust, closures fall into three kinds: FnOnce, FnMut and Fn.
// FnOnce assumed to be able to be called just once if it is not Clone. It is impossible to
// write recursive FnOnce that is not Clone.
// All FnMut are also FnOnce, although you can call them multiple times, they are not allow to
// have a reference to themselves. So it is also not possible to write recursive FnMut closures
// that is not Clone.
// All Fn are also FnMut, and all closures of Fn are also Clone. However, programmers can create
// Fn objects that are not Clone
 
// This will work for all Fn objects, not just closures
// And it is a little bit more efficient for Fn closures as it do not clone itself.
impl<T, R, F> Apply<T, R> for F where F:
  Fn(&dyn Apply<T, R>, T) -> R
{
  fn apply(&self, f: &dyn Apply<T, R>, t: T) -> R {
    self(f, t)
 
    // NOTE: Each letter is an individual symbol.
    // (λx.(λy.xxy))(λx.(λy.f(λz.xxz)y))t
    // => (λx.xx)(λx.f(xx))t
    // => (Yf)t
  }
}
 
// This works for all closures that is Clone, and those are Fn.
// impl<T, R, F> Apply<T, R> for F where F: FnOnce( &Apply<T, R>, T ) -> R + Clone {
//     fn apply( &self, f: &Apply<T, R>, t: T ) -> R {
//         (self.clone())( f, t )
 
//         // If we were to pass in self as f, we get -
//         // NOTE: Each letter is an individual symbol.
//         // λf.λt.sft
//         // => λs.λt.sst [s/f]
//         // => λs.ss
//     }
// }
 
// Before 1.26 we have some limitations and so we need some workarounds. But now impl Trait is stable and we can
// write the following:
 
fn y<T,R>(f:impl Fn(&dyn Fn(T) -> R, T) -> R) -> impl Fn(T) -> R {
  move |t| (
    |x: &dyn Apply<T,R>, y| x.apply(x, y)
  ) (
    &|x: &dyn Apply<T,R>, y| f(
      &|z| x.apply(x,z),
      y
    ),
    t
  )
}
 
// fn y<T,R>(f:impl FnOnce(&Fn(T) -> R, T) -> R + Clone) -> impl FnOnce(T) -> R {
//    |t| (|x: &Apply<T,R>,y| x.apply(x,y))
//        (&move |x:&Apply<T,R>,y| f(&|z| x.apply(x,z), y), t)
 
//     // NOTE: Each letter is an individual symbol.
//     // (λx.(λy.xxy))(λx.(λy.f(λz.xxz)y))t
//     // => (λx.xx)(λx.f(xx))t
//     // => (Yf)t
// }
 
// Previous version removed as they are just hacks when impl Trait is not available.
 
fn fac(n: usize) -> usize {
  let almost_fac = |f: &dyn Fn(usize) -> usize, x|
    if x == 0 {
      1
    } else {
      x * f(x - 1)
    }
  ;
  let fac = y( almost_fac );
  fac(n)
}
 
fn fib( n: usize ) -> usize {
  let almost_fib = |f: &dyn Fn(usize) -> usize, x|
    if x < 2 {
      1
    } else {
      f(x - 2) + f(x - 1)
    };
  let fib = y(almost_fib);
  fib(n)
}
 
fn optimal_fib( n: usize ) -> usize {
  let almost_fib = |f: &dyn Fn((usize,usize,usize)) -> usize, (i0,i1,x)| 
    match x {
      0 => i0,
      1 => i1,
      x => f((i1,i0+i1, x-1))
    }        
  ;
  let fib = |x| y(almost_fib)((1,1,x));
  fib(n)
}
 
fn test_y() {
  println!("{}", fac(10));
  println!("{}", fib(10));
  println!("{}", optimal_fib(10));
}
